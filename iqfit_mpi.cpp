// iqfit_mpi.cpp (Translated to English, with meaningful variable names and enhanced comments)
// MPI-based parallel solver for the IQ-Fit puzzle (11x5 board, 12 unique pieces).
// Each MPI rank explores a disjoint set of possible placements for the first piece.

#include <mpi.h>
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <numeric>
#include <array>

// Board and puzzle parameters
constexpr int BOARD_WIDTH = 11;
constexpr int BOARD_HEIGHT = 5;
constexpr int TOTAL_CELLS = BOARD_WIDTH * BOARD_HEIGHT;
constexpr int TOTAL_PIECES = 12;

// Each shape string defines a base piece using "xy" format 
static const std::vector<std::string> basePieceShapes = {
    "01 10 11 21 31", "01 10 11 21 22", "10 11 12 13 03",
    "01 11 10 02", "00 01 02 12 13", "02 12 11 21 20",
    "02 12 11 10", "02 12 22 21 20", "01 11 10",
    "01 02 11 12 10", "01 11 10 21", "00 01 11 21 20"
};

// Precomputed bitmask placements for each piece and its variations
std::vector<std::vector<uint64_t>> piecePlacementMasks(TOTAL_PIECES);
// List of board cell indices covered by each valid placement
std::vector<std::vector<std::vector<int>>> piecePlacementCells(TOTAL_PIECES);
// For each piece and board cell: which placements cover that cell
std::vector<std::vector<std::vector<int>>> piecePlacementsByCell(TOTAL_PIECES, std::vector<std::vector<int>>(TOTAL_CELLS));

// Representation of the board as a 1D character array
using BoardRepresentation = std::array<char, TOTAL_CELLS>;

// Parse a piece shape string into a list of coordinate pairs
static std::vector<std::pair<int,int>> parsePieceShape(const std::string &shapeStr) {
    std::vector<std::pair<int,int>> coordinates;
    for (size_t i = 0; i + 1 < shapeStr.size(); i += 3) {
        int x = shapeStr[i] - '0';
        int y = shapeStr[i+1] - '0';
        coordinates.emplace_back(x, y);
    }
    return coordinates;
}

// Generate all unique orientations (rotations + reflections) of a piece
static std::vector<std::vector<std::pair<int,int>>> generateUniqueOrientations(const std::vector<std::pair<int,int>> &baseCoords) {
    std::set<std::vector<std::pair<int,int>>> uniqueOrientations;
    for (int reflect = 0; reflect < 2; ++reflect) {
        for (int rot = 0; rot < 4; ++rot) {
            std::vector<std::pair<int,int>> transformed;
            for (const auto &coord : baseCoords) {
                int x = reflect ? -coord.first : coord.first;
                int y = coord.second;
                for (int r = 0; r < rot; ++r) {
                    int temp = x;
                    x = y;
                    y = -temp;
                }
                transformed.emplace_back(x, y);
            }
            // Normalize to top-left origin
            int minX = INT32_MAX, minY = INT32_MAX;
            for (const auto &p : transformed) {
                minX = std::min(minX, p.first);
                minY = std::min(minY, p.second);
            }
            for (auto &p : transformed) {
                p.first -= minX;
                p.second -= minY;
            }
            std::sort(transformed.begin(), transformed.end());
            uniqueOrientations.insert(transformed);
        }
    }
    return std::vector<std::vector<std::pair<int,int>>>(uniqueOrientations.begin(), uniqueOrientations.end());
}

// Precompute all legal placements for every piece in all orientations
static void precomputeAllPiecePlacements() {
    for (int pieceIdx = 0; pieceIdx < TOTAL_PIECES; ++pieceIdx) {
        auto baseCoords = parsePieceShape(basePieceShapes[pieceIdx]);
        auto allOrientations = generateUniqueOrientations(baseCoords);

        for (const auto &shape : allOrientations) {
            int maxX = 0, maxY = 0;
            for (const auto &coord : shape) {
                maxX = std::max(maxX, coord.first);
                maxY = std::max(maxY, coord.second);
            }
            int shapeWidth = maxX + 1;
            int shapeHeight = maxY + 1;

            for (int yOffset = 0; yOffset <= BOARD_HEIGHT - shapeHeight; ++yOffset) {
                for (int xOffset = 0; xOffset <= BOARD_WIDTH - shapeWidth; ++xOffset) {
                    uint64_t placementMask = 0ULL;
                    std::vector<int> cellIndices;
                    bool validPlacement = true;
                    for (const auto &coord : shape) {
                        int x = xOffset + coord.first;
                        int y = yOffset + coord.second;
                        int cellIdx = y * BOARD_WIDTH + x;
                        if (cellIdx < 0 || cellIdx >= TOTAL_CELLS) {
                            validPlacement = false;
                            break;
                        }
                        placementMask |= (1ULL << cellIdx);
                        cellIndices.push_back(cellIdx);
                    }
                    if (!validPlacement) continue;
                    int placementIdx = piecePlacementMasks[pieceIdx].size();
                    piecePlacementMasks[pieceIdx].push_back(placementMask);
                    piecePlacementCells[pieceIdx].push_back(cellIndices);
                    for (int cell : cellIndices) {
                        piecePlacementsByCell[pieceIdx][cell].push_back(placementIdx);
                    }
                }
            }
        }
    }
}

// Recursive backtracking search to find valid solutions
static void recursiveSolver(
    uint64_t currentBoardMask,
    std::array<bool, TOTAL_PIECES> &usedPieces,
    BoardRepresentation &currentBoard,
    std::vector<BoardRepresentation> &foundSolutions
) {
    // Base case: all pieces placed
    if (std::all_of(usedPieces.begin(), usedPieces.end(), [](bool used) { return used; })) {
        foundSolutions.push_back(currentBoard);
        return;
    }

    // Find the first empty cell
    int firstEmptyCell = 0;
    while (firstEmptyCell < TOTAL_CELLS && ((currentBoardMask >> firstEmptyCell) & 1ULL)) {
        ++firstEmptyCell;
    }
    if (firstEmptyCell >= TOTAL_CELLS) return;

    // Try all unused pieces that can cover the current cell
    for (int pieceIdx = 0; pieceIdx < TOTAL_PIECES; ++pieceIdx) {
        if (usedPieces[pieceIdx]) continue;
        for (int placementIdx : piecePlacementsByCell[pieceIdx][firstEmptyCell]) {
            uint64_t placementMask = piecePlacementMasks[pieceIdx][placementIdx];
            if ((placementMask & currentBoardMask) != 0ULL) continue;

            // Place the piece
            usedPieces[pieceIdx] = true;
            uint64_t newMask = currentBoardMask | placementMask;
            for (int cell : piecePlacementCells[pieceIdx][placementIdx]) {
                currentBoard[cell] = char('A' + pieceIdx);
            }
            recursiveSolver(newMask, usedPieces, currentBoard, foundSolutions);
            // Backtrack
            usedPieces[pieceIdx] = false;
            for (int cell : piecePlacementCells[pieceIdx][placementIdx]) {
                currentBoard[cell] = '.';
            }
        }
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int totalRanks, rankId;
    MPI_Comm_size(MPI_COMM_WORLD, &totalRanks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankId);

    double startTime = MPI_Wtime();
    precomputeAllPiecePlacements();

    int totalStartingPlacements = piecePlacementMasks[0].size();
    std::vector<BoardRepresentation> localSolutions;
    BoardRepresentation initialBoard;
    initialBoard.fill('.');
    std::array<bool, TOTAL_PIECES> initialUsed;
    initialUsed.fill(false);

    // Distribute first-piece placements among MPI ranks
    for (int i = rankId; i < totalStartingPlacements; i += totalRanks) {
        BoardRepresentation currentBoard = initialBoard;
        auto used = initialUsed;
        uint64_t currentMask = piecePlacementMasks[0][i];
        used[0] = true;
        for (int cell : piecePlacementCells[0][i]) {
            currentBoard[cell] = 'A';
        }
        recursiveSolver(currentMask, used, currentBoard, localSolutions);
    }

    // Collect solution counts
    int localCount = localSolutions.size();
    std::vector<int> solutionCounts;
    if (rankId == 0) {
        solutionCounts.resize(totalRanks);
    }
    MPI_Gather(&localCount, 1, MPI_INT,
               solutionCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Flatten local solutions to char buffer
    int localChars = localCount * TOTAL_CELLS;
    std::vector<char> localBuffer(localChars);
    for (int i = 0; i < localCount; ++i) {
        std::memcpy(&localBuffer[i * TOTAL_CELLS], localSolutions[i].data(), TOTAL_CELLS);
    }

    // Setup receive buffers on rank 0
    std::vector<int> recvCounts, displacements;
    std::vector<char> allSolutionsBuffer;
    if (rankId == 0) {
        recvCounts.resize(totalRanks);
        displacements.resize(totalRanks);
        int offset = 0;
        for (int i = 0; i < totalRanks; ++i) {
            recvCounts[i] = solutionCounts[i] * TOTAL_CELLS;
            displacements[i] = offset;
            offset += recvCounts[i];
        }
        allSolutionsBuffer.resize(offset);
    }

    // Gather all boards into rank 0
    MPI_Gatherv(localBuffer.data(), localChars, MPI_CHAR,
                allSolutionsBuffer.data(), recvCounts.data(), displacements.data(),
                MPI_CHAR, 0, MPI_COMM_WORLD);

    // Output results to file from rank 0
    if (rankId == 0) {
        std::ofstream outputFile("solutions.txt");
        if (!outputFile.is_open()) {
            std::cerr << "Error: Could not open solutions.txt\n";
        } else {
            int totalSolutions = std::accumulate(solutionCounts.begin(), solutionCounts.end(), 0);
            for (int r = 0; r < totalRanks; ++r) {
                int count = solutionCounts[r];
                for (int s = 0; s < count; ++s) {
                    const char *boardData = allSolutionsBuffer.data() + displacements[r] + s * TOTAL_CELLS;
                    for (int row = 0; row < BOARD_HEIGHT; ++row) {
                        outputFile.write(boardData + row * BOARD_WIDTH, BOARD_WIDTH);
                        outputFile.put('\n');
                    }
                    outputFile.put('\n');
                }
            }
            outputFile.close();
            std::cout << "Total solutions: " << totalSolutions << "\n";
        }
        double endTime = MPI_Wtime();
        std::cout << "Elapsed time: " << (endTime - startTime) << " seconds\n";
    }

    MPI_Finalize();
    return 0;
}
