# Shamir's Secret Sharing Solver

A high-performance Java implementation of Shamir's Secret Sharing scheme with intelligent algorithm selection and comprehensive performance analysis.

[![Java](https://img.shields.io/badge/Java-11+-orange.svg)](https://www.oracle.com/java/)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Performance](https://img.shields.io/badge/Performance-Optimized-green.svg)]()

## 🚀 Features

- **Dual Algorithm Implementation**: Lagrange Interpolation & Gaussian Elimination
- **Intelligent Algorithm Selection**: Automatic optimization based on problem size
- **Modular Arithmetic**: Secure finite field operations using prime `2^127 - 1`
- **Performance Benchmarking**: Comprehensive analysis tools with scaling studies
- **Cross-Validation**: Multiple point combinations for result verification
- **Production Ready**: Extensive error handling and numerical stability testing
- **Parallel Processing**: Multi-threaded support for large computations

## 📊 Performance Highlights

| Problem Size (k) | Lagrange Time | Gaussian Time | Speed Advantage |
|------------------|---------------|---------------|-----------------|
| k = 10          | 0.093 ms      | 0.226 ms      | **2.4x faster** |
| k = 20          | 0.236 ms      | 1.260 ms      | **5.3x faster** |
| k = 30          | 0.391 ms      | 4.048 ms      | **10.4x faster** |

## 🏗️ Architecture

### Three Implementation Levels:

1. **`ShamirSecretSolver.java`** - Basic Lagrange-only implementation
2. **`GaussianShamirSecretSolver.java`** - Production-grade hybrid solver ⭐
3. **`PerformanceTester.java`** - Comprehensive benchmarking suite

## 🚦 Quick Start

### Prerequisites

- Java 11 or higher
- Jackson JSON library for JSON parsing
- Maven or Gradle (for dependency management)

### Dependencies

```xml
<dependency>
    <groupId>com.fasterxml.jackson.core</groupId>
    <artifactId>jackson-databind</artifactId>
    <version>2.15.2</version>
</dependency>
```

### Basic Usage

```java
// Initialize the solver
GaussianShamirSecretSolver solver = new GaussianShamirSecretSolver();

// Solve from JSON file
BigInteger secret = solver.solveFromFile("testcase1.json");
System.out.println("Recovered Secret: " + secret);
```

### Input Format

```json
{
  "keys": {
    "n": 4,
    "k": 3
  },
  "1": {
    "base": "10",
    "value": "4"
  },
  "2": {
    "base": "2", 
    "value": "111"
  },
  "3": {
    "base": "10",
    "value": "12"
  }
}
```

## 🎯 Algorithm Selection Strategy

The hybrid solver automatically chooses the optimal algorithm:

```java
if (k <= 10) {
    // Use Lagrange Interpolation - O(k²)
    return lagrangeInterpolation(points);
} else {
    // Use Gaussian Elimination - O(k³) 
    return gaussianEliminationMethod(points);
}
```

### When to Use Each Algorithm:

| Algorithm | Best For | Complexity | Memory |
|-----------|----------|------------|---------|
| **Lagrange** | k ≤ 20, Real-time apps | O(k²) | O(k) |
| **Gaussian** | k > 20, Matrix operations | O(k³) | O(k²) |
| **Hybrid** | Production systems | Optimal | Adaptive |

## 🔬 Performance Analysis

### Run Comprehensive Benchmarks

```bash
javac PerformanceTester.java
java PerformanceTester
```

### Sample Output:

```
=== SHAMIR'S SECRET SHARING PERFORMANCE BENCHMARK ===
| k Size | Lagrange Time | Gaussian Time | Memory (MB) | Recommended |
|--------|---------------|---------------|-------------|-------------|
|     10 |       0.093 ms |       0.226 ms |      0.00 |    Lagrange |
|     20 |       0.236 ms |       1.260 ms |      0.01 |    Lagrange |
|     30 |       0.391 ms |       4.048 ms |      0.03 |    Lagrange |

=== SCALING ANALYSIS ===
Lagrange Scaling (should be ~O(k²)):
k=10 to k=20: actual ratio=2.54, expected O(k²)=4.00
k=20 to k=40: actual ratio=1.66, expected O(k²)=4.00

=== NUMERICAL STABILITY TEST ===
Test Case              | Lagrange Error | Gaussian Error | Winner
-----------------------|----------------|----------------|--------
Sequential points      |        Perfect |        Perfect | Lagrange
Large gaps             |        Perfect |        Perfect | Lagrange
Random distribution    |        Perfect |        Perfect | Lagrange
```

## 🛡️ Security Features

- **Finite Field Arithmetic**: All operations in GF(p) using Mersenne prime
- **Modular Inverse**: Extended Euclidean algorithm implementation  
- **Overflow Protection**: Safe BigInteger operations
- **Input Validation**: Comprehensive base and format checking

## 📈 Use Cases

### Cryptographic Applications
- Distributed key management
- Multi-party computation protocols
- Secure backup systems
- Blockchain wallet recovery

### Performance-Critical Systems
- Real-time secret reconstruction
- High-throughput cryptographic services
- Embedded security systems
- IoT device authentication

## 🧪 Testing

### Run Basic Tests
```bash
# Compile all classes
javac -cp ".:jackson-databind-2.15.2.jar:jackson-core-2.15.2.jar:jackson-annotations-2.15.2.jar" *.java

# Run main solver
java -cp ".:jackson-databind-2.15.2.jar:jackson-core-2.15.2.jar:jackson-annotations-2.15.2.jar" GaussianShamirSecretSolver

# Run performance benchmarks
java PerformanceTester
```

### Test Cases Included
- Small polynomials (k=3-5)
- Medium complexity (k=10-20) 
- Large systems (k=25-30)
- Edge cases (sequential points, large gaps, small differences)

## 🔧 Configuration

### Adjustable Parameters

```java
// Field prime for modular arithmetic
private static final BigInteger FIELD_PRIME = 
    new BigInteger("170141183460469231731687303715884105727");

// Algorithm selection threshold
private static final int GAUSSIAN_THRESHOLD = 10;

// Performance testing iterations
private static final int TEST_ITERATIONS = 1000;
```

## 📊 Complexity Analysis

| Operation | Lagrange | Gaussian | Memory |
|-----------|----------|----------|---------|
| **Time** | O(k²) | O(k³) | - |
| **Space** | O(k) | O(k²) | BigInteger overhead |
| **Stability** | Excellent | Excellent | IEEE 754 precision |

## 🤝 Contributing

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Run performance tests (`java PerformanceTester`)
4. Commit your changes (`git commit -m 'Add amazing feature'`)
5. Push to the branch (`git push origin feature/amazing-feature`)
6. Open a Pull Request

### Development Guidelines
- Maintain O(k²) performance for small k
- Add comprehensive test cases
- Update benchmarks for new features
- Follow existing code style

## 📝 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🏆 Acknowledgments

- **Adi Shamir** - Original Secret Sharing scheme (1979)
- **Joseph-Louis Lagrange** - Polynomial interpolation method
- **Carl Friedrich Gauss** - Gaussian elimination algorithm
- **Performance Optimization** - Modern JVM hotspot optimizations

## 📚 References

- [Shamir's Secret Sharing - Original Paper](https://dl.acm.org/doi/10.1145/359168.359176)
- [Lagrange Interpolation](https://en.wikipedia.org/wiki/Lagrange_polynomial)
- [Gaussian Elimination](https://en.wikipedia.org/wiki/Gaussian_elimination)
- [Finite Field Arithmetic](https://en.wikipedia.org/wiki/Finite_field_arithmetic)

## 📞 Support

- **Issues**: [GitHub Issues](https://github.com/yourusername/shamir-secret-sharing/issues)
- **Discussions**: [GitHub Discussions](https://github.com/yourusername/shamir-secret-sharing/discussions)
- **Performance Questions**: Check the benchmark results in `PerformanceTester.java`

---

⭐ **Star this repository if you find it useful!** ⭐
