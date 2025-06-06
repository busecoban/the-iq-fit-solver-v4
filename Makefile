# Makefile for IQ Fit MPI Solver - Multi-core Test Ready

CXX = mpic++
CXXFLAGS = -std=c++11 -O3
TARGET = iqfit_mpi
SRC = iqfit_mpi.cpp

# Build target
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

# Run targets with different core counts
run1: $(TARGET)
	@echo "🚀 Running with 1 core..."
	@mkdir -p log
	mpirun -np 1 ./$(TARGET) | tee log/run1.txt

run2: $(TARGET)
	@echo "🚀 Running with 2 core..."
	@mkdir -p log
	mpirun -np 2 ./$(TARGET) | tee log/run2.txt

run4: $(TARGET)
	@echo "🚀 Running with 4 cores..."
	@mkdir -p log
	mpirun -np 4 ./$(TARGET) | tee log/run4.txt

run8: $(TARGET)
	@echo "🚀 Running with 8 cores..."
	@mkdir -p log
	mpirun -np 8 ./$(TARGET) | tee log/run8.txt

run12: $(TARGET)
	@echo "🚀 Running with 12 cores..."
	@mkdir -p log
	mpirun -np 12 ./$(TARGET) | tee log/run12.txt

# Clean build and output files
clean:
	rm -f $(TARGET) solutions.txt
	rm -rf log
