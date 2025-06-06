// iqfit_mpi.cpp
// MPI-based parallel solver for the IQ-Fit puzzle (11×5 board, 12 pieces).
// C++11 ile uyumlu hale getirildi, ayrıca çalışma süresi MPI_Wtime() ile ölçülür.

#include <mpi.h>
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <numeric>   // std::accumulate için
#include <array>     // std::array kullanımı için

constexpr int BOARD_WIDTH = 11;
constexpr int BOARD_HEIGHT = 5;
constexpr int NUM_CELLS = BOARD_WIDTH * BOARD_HEIGHT; // 55
constexpr int NUM_PIECES = 12;

// Ham parça tanımları: her string "xy" biçiminde koordinatlardan oluşur
static const std::vector<std::string> pieceShapes = {
    "01 10 11 21 31",
    "01 10 11 21 22",
    "10 11 12 13 03",
    "01 11 10 02",
    "00 01 02 12 13",
    "02 12 11 21 20",
    "02 12 11 10",
    "02 12 22 21 20",
    "01 11 10",
    "01 02 11 12 10",
    "01 11 10 21",
    "00 01 11 21 20"
};

// Her parça p için: yerleşim maskeleri ve ilgili hücre indeksleri
std::vector<std::vector<uint64_t>> placementMasks(NUM_PIECES);
std::vector<std::vector<std::vector<int>>> placementCells(NUM_PIECES);

// Her parça p ve her hücre c için: c hücresini kapsayan yerleşim-indeksleri
std::vector<std::vector<std::vector<int>>> placementsByCell(NUM_PIECES, std::vector<std::vector<int>>(NUM_CELLS));

// Bir tahtanın karakter tamponu (55 hücre)
using BoardChars = std::array<char, NUM_CELLS>;

// Bir shape string'ini (örneğin "01 10 11 21 31") listeye çevirir
static std::vector<std::pair<int,int>> parseShape(const std::string &s) {
    std::vector<std::pair<int,int>> coords;
    for (size_t i = 0; i + 1 < s.size(); i += 3) {
        int x = s[i] - '0';
        int y = s[i+1] - '0';
        coords.emplace_back(x, y);
    }
    return coords;
}

// Bir base koordinat dizisinden (örneğin {(0,1),(1,0),...})
// tüm benzersiz oryantasyonları üretir (rotasyon + yansıma)
static std::vector<std::vector<std::pair<int,int>>> generateOrientations(const std::vector<std::pair<int,int>> &base) {
    std::set<std::vector<std::pair<int,int>>> uniqueShapes;

    for (int reflect = 0; reflect < 2; ++reflect) {
        for (int rot = 0; rot < 4; ++rot) {
            std::vector<std::pair<int,int>> transformed;
            // Önce yansıma (x ekseninde) uygula
            for (const auto &p0 : base) {
                int x0 = p0.first;
                int y0 = p0.second;
                int x = (reflect ? -x0 : x0);
                int y = y0;
                // Rotasyon rot kez 90° döndür
                for (int r = 0; r < rot; ++r) {
                    int tx = x;
                    x = y;
                    y = -tx;
                }
                transformed.emplace_back(x, y);
            }
            // Normalize et: min x,y = 0 olacak şekilde kaydır
            int minX = INT32_MAX, minY = INT32_MAX;
            for (size_t i = 0; i < transformed.size(); ++i) {
                minX = std::min(minX, transformed[i].first);
                minY = std::min(minY, transformed[i].second);
            }
            for (size_t i = 0; i < transformed.size(); ++i) {
                transformed[i].first  -= minX;
                transformed[i].second -= minY;
            }
            // Sıralayarak benzersizleştir
            std::sort(transformed.begin(), transformed.end());
            uniqueShapes.insert(transformed);
        }
    }

    std::vector<std::vector<std::pair<int,int>>> result;
    for (auto &shape : uniqueShapes) {
        result.push_back(shape);
    }
    return result;
}

// Her parça için tüm geçerli yerleşimleri önceden hesapla
static void precomputePlacements() {
    for (int p = 0; p < NUM_PIECES; ++p) {
        // Base shape'i parse et
        auto baseCoords = parseShape(pieceShapes[p]);
        // Bütün benzersiz oryantasyonları al
        auto orientations = generateOrientations(baseCoords);

        // Her oryantasyon için tüm kaydırmaları dene
        for (const auto &shape : orientations) {
            // Oriented shape'in genişlik ve yüksekliğini bul
            int maxX = 0, maxY = 0;
            for (size_t i = 0; i < shape.size(); ++i) {
                int x = shape[i].first;
                int y = shape[i].second;
                maxX = std::max(maxX, x);
                maxY = std::max(maxY, y);
            }
            int shapeW = maxX + 1;
            int shapeH = maxY + 1;

            // Tahtanın içinde kaydırarak yerleştir
            for (int oy = 0; oy <= BOARD_HEIGHT - shapeH; ++oy) {
                for (int ox = 0; ox <= BOARD_WIDTH  - shapeW; ++ox) {
                    uint64_t mask = 0ULL;
                    std::vector<int> cells;
                    bool valid = true;
                    for (size_t i = 0; i < shape.size(); ++i) {
                        int x = ox + shape[i].first;
                        int y = oy + shape[i].second;
                        int idx = y * BOARD_WIDTH + x;
                        if (idx < 0 || idx >= NUM_CELLS) {
                            valid = false;
                            break;
                        }
                        mask |= (1ULL << idx);
                        cells.push_back(idx);
                    }
                    if (!valid) continue;
                    // Geçerli yerleşimi kaydet
                    int placementIndex = placementMasks[p].size();
                    placementMasks[p].push_back(mask);
                    placementCells[p].push_back(cells);
                    // placementsByCell[p][c] vektörüne ekle
                    for (size_t j = 0; j < cells.size(); ++j) {
                        int idx = cells[j];
                        placementsByCell[p][idx].push_back(placementIndex);
                    }
                }
            }
        }
    }
}

// Düğüme dayalı geri izleme (backtracking):
// Mevcut boardMask, hangi parçaların kullanıldığı (used[]), ve boardChars tamponu
// çözümleri solutions vektörüne ekler.
static void search(
    uint64_t boardMask,
    std::array<bool, NUM_PIECES> &used,
    BoardChars &boardChars,
    std::vector<BoardChars> &solutions
) {
    // Tüm parçalar kullanıldı mı?
    bool done = true;
    for (int i = 0; i < NUM_PIECES; ++i) {
        if (!used[i]) {
            done = false;
            break;
        }
    }
    if (done) {
        solutions.push_back(boardChars);
        return;
    }

    // İlk boş hücreyi bul (bit = 0)
    int c = 0;
    while (c < NUM_CELLS && ((boardMask >> c) & 1ULL)) {
        ++c;
    }
    if (c >= NUM_CELLS) return; // Tüm hücreler dolu ama parçalar bitmedi => çık

    // c hücresini kapsayan her kullanılmamış parçayı dene
    for (int p = 0; p < NUM_PIECES; ++p) {
        if (used[p]) continue;
        // Bu parça için c hücresini kapsayan tüm yerleşim-indeksleri
        for (int idx : placementsByCell[p][c]) {
            uint64_t pmask = placementMasks[p][idx];
            // Örtüşme var mı?
            if ((pmask & boardMask) != 0ULL) continue;
            // Parçayı yerleştir
            used[p] = true;
            uint64_t newMask = boardMask | pmask;
            for (size_t j = 0; j < placementCells[p][idx].size(); ++j) {
                int cell = placementCells[p][idx][j];
                boardChars[cell] = char('A' + p);
            }
            // Rekürsif çağrı
            search(newMask, used, boardChars, solutions);
            // Geri al (backtrack)
            used[p] = false;
            for (size_t j = 0; j < placementCells[p][idx].size(); ++j) {
                int cell = placementCells[p][idx][j];
                boardChars[cell] = '.';
            }
        }
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int world_size = 1, world_rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Zaman ölçümü başlat
    double startTime = MPI_Wtime();

    // Tüm parçalar için yerleşimleri önceden hesapla (tüm ranklerde)
    precomputePlacements();

    // Parça 0 (A) yerleşim sayısı
    int totalPlacements0 = static_cast<int>(placementMasks[0].size());

    // Bu rank'in bulduğu çözümler
    std::vector<BoardChars> localSolutions;

    // Boş tahta karakter tamponu ('.' ile doldurulmuş)
    BoardChars boardCharsInit;
    boardCharsInit.fill('.');

    // Hangi parçalar kullanıldı?
    std::array<bool, NUM_PIECES> usedInit;
    usedInit.fill(false);

    // Her rank, parça 0'ın farklı indekslerini işlesin: i0 = rank, rank + world_size, ...
    for (int i0 = world_rank; i0 < totalPlacements0; i0 += world_size) {
        // Her iterasyonda tahtayı temizle
        BoardChars boardChars = boardCharsInit;
        auto used = usedInit;
        uint64_t boardMask = 0ULL;

        // Parça 0'ı i0 konumunda yerleştir
        used[0] = true;
        boardMask = placementMasks[0][i0];
        for (size_t j = 0; j < placementCells[0][i0].size(); ++j) {
            int cell = placementCells[0][i0][j];
            boardChars[cell] = 'A';
        }

        // Kalan parçaları arama
        search(boardMask, used, boardChars, localSolutions);
        // Bir sonraki i0 için aynı başlangıç adımları yeniden oluşturulacak
    }

    // Şimdi local çözümleri rank 0'a toplayalım
    int localCount = static_cast<int>(localSolutions.size());
    std::vector<int> allCounts;
    if (world_rank == 0) {
        allCounts.resize(world_size, 0);
    }
    MPI_Gather(&localCount, 1, MPI_INT,
               (world_rank == 0 ? allCounts.data() : nullptr),
               1, MPI_INT, 0, MPI_COMM_WORLD);

    // localSolutions[i] her biri 55 karakterlik, tüm rank'ların flatten edilmiş hali
    int localChars = localCount * NUM_CELLS;
    std::vector<char> localBuf(localChars);
    for (int i = 0; i < localCount; ++i) {
        std::memcpy(&localBuf[i * NUM_CELLS], localSolutions[i].data(), NUM_CELLS);
    }

    // Rank 0'de alıcı tamponları ayarla
    std::vector<int> recvCounts;
    std::vector<int> displs;
    std::vector<char> recvBuf;
    if (world_rank == 0) {
        recvCounts.resize(world_size);
        displs.resize(world_size);
        int offset = 0;
        for (int i = 0; i < world_size; ++i) {
            recvCounts[i] = allCounts[i] * NUM_CELLS;
            displs[i]    = offset;
            offset       += recvCounts[i];
        }
        recvBuf.resize(offset);
    }

    // Her rank'in buffer'ını rank 0'a topla
    MPI_Gatherv(localBuf.data(), localChars, MPI_CHAR,
                (world_rank == 0 ? recvBuf.data()    : nullptr),
                (world_rank == 0 ? recvCounts.data() : nullptr),
                (world_rank == 0 ? displs.data()     : nullptr),
                MPI_CHAR, 0, MPI_COMM_WORLD);

    // Rank 0 tüm çözümleri "solutions.txt" dosyasına yazsın
    if (world_rank == 0) {
        std::ofstream ofs("solutions.txt");
        if (!ofs.is_open()) {
            std::cerr << "Hata: solutions.txt açılamadı\n";
        } else {
            int totalSolutions = std::accumulate(allCounts.begin(), allCounts.end(), 0);
            for (int r = 0; r < world_size; ++r) {
                int count = allCounts[r];
                for (int s = 0; s < count; ++s) {
                    // Bir çözüme ait 55 karaktere işaret eden pointer
                    const char *boardPtr = recvBuf.data() + displs[r] + s * NUM_CELLS;
                    // 5 satır × 11 sütun yaz
                    for (int row = 0; row < BOARD_HEIGHT; ++row) {
                        ofs.write(boardPtr + row * BOARD_WIDTH, BOARD_WIDTH);
                        ofs.put('\n');
                    }
                    ofs.put('\n');
                }
            }
            ofs.close();
            std::cout << "Toplam çözüm sayısı: " << totalSolutions << "\n";
        }

        // Zaman ölçümünü al
        double endTime = MPI_Wtime();
        double elapsed = endTime - startTime;
        std::cout << "Toplam çalışma süresi: " << elapsed << " saniye\n";
    }

    MPI_Finalize();
    return 0;
} 