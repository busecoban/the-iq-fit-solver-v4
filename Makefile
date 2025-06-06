CXX = mpic++
CXXFLAGS = -std=c++11 -O3
TARGET = iqfit_mpi
SRC = iqfit_mpi.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

run: $(TARGET)
	@echo "🚀 12 çekirdek ile çalıştırılıyor..."
	mpirun -np 12 ./$(TARGET)

clean:
	rm -f $(TARGET) solutions.txt
