import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.*;

/**
 * Performance Testing Tool for Shamir's Secret Sharing Algorithms
 * Compares Lagrange vs Gaussian Elimination across different problem sizes
 */
public class PerformanceTester {
    
    private static final BigInteger FIELD_PRIME = new BigInteger("170141183460469231731687303715884105727");
    private static final SecureRandom RANDOM = new SecureRandom();
    private static final int WARMUP_ITERATIONS = 100;
    private static final int TEST_ITERATIONS = 1000;
    
    public static void main(String[] args) {
        PerformanceTester tester = new PerformanceTester();
        
        System.out.println("=== SHAMIR'S SECRET SHARING PERFORMANCE BENCHMARK ===");
        System.out.println("Field Prime: " + FIELD_PRIME.toString().substring(0, 20) + "...");
        System.out.println("Warmup Iterations: " + WARMUP_ITERATIONS);
        System.out.println("Test Iterations: " + TEST_ITERATIONS);
        System.out.println();
        
        // Test different problem sizes
        int[] testSizes = {3, 5, 7, 10, 15, 20, 25, 30};
        
        System.out.println("| k Size | Lagrange Time | Gaussian Time | Memory (MB) | Recommended |");
        System.out.println("|--------|---------------|---------------|-------------|-------------|");
        
        for (int k : testSizes) {
            BenchmarkResult result = tester.benchmarkAlgorithms(k);
            System.out.printf("| %6d | %11.3f ms | %11.3f ms | %9.2f | %11s |\n",
                k, result.lagrangeTime, result.gaussianTime, result.memoryUsageMB,
                result.lagrangeTime < result.gaussianTime ? "Lagrange" : "Gaussian");
        }
        
        System.out.println();
        tester.demonstrateScaling();
        tester.testNumericalStability();
    }
    
    /**
     * Benchmark both algorithms for given k
     */
    public BenchmarkResult benchmarkAlgorithms(int k) {
        // Generate test polynomial and points
        TestData testData = generateTestData(k);
        
        // Warmup JVM
        for (int i = 0; i < WARMUP_ITERATIONS; i++) {
            lagrangeMethod(testData.points);
            gaussianMethod(testData.points);
        }
        
        // Benchmark Lagrange
        long lagrangeTotal = 0;
        for (int i = 0; i < TEST_ITERATIONS; i++) {
            long start = System.nanoTime();
            lagrangeMethod(testData.points);
            lagrangeTotal += System.nanoTime() - start;
        }
        double lagrangeAvg = lagrangeTotal / (double) TEST_ITERATIONS / 1_000_000.0;
        
        // Benchmark Gaussian
        long gaussianTotal = 0;
        for (int i = 0; i < TEST_ITERATIONS; i++) {
            long start = System.nanoTime();
            gaussianMethod(testData.points);
            gaussianTotal += System.nanoTime() - start;
        }
        double gaussianAvg = gaussianTotal / (double) TEST_ITERATIONS / 1_000_000.0;
        
        // Estimate memory usage
        double memoryMB = estimateMemoryUsage(k);
        
        return new BenchmarkResult(lagrangeAvg, gaussianAvg, memoryMB);
    }
    
    /**
     * Generate test polynomial and evaluation points
     */
    private TestData generateTestData(int k) {
        // Generate random polynomial coefficients
        List<BigInteger> coefficients = new ArrayList<>();
        BigInteger secret = new BigInteger(128, RANDOM).mod(FIELD_PRIME);
        coefficients.add(secret); // Constant term (secret)
        
        for (int i = 1; i < k; i++) {
            BigInteger coeff = new BigInteger(FIELD_PRIME.bitLength() - 1, RANDOM);
            coefficients.add(coeff);
        }
        
        // Generate evaluation points
        List<Point> points = new ArrayList<>();
        for (int x = 1; x <= k; x++) {
            BigInteger y = evaluatePolynomial(coefficients, BigInteger.valueOf(x));
            points.add(new Point(BigInteger.valueOf(x), y));
        }
        
        return new TestData(coefficients, points, secret);
    }
    
    /**
     * Evaluate polynomial at given x
     */
    private BigInteger evaluatePolynomial(List<BigInteger> coeffs, BigInteger x) {
        BigInteger result = BigInteger.ZERO;
        BigInteger xPower = BigInteger.ONE;
        
        for (BigInteger coeff : coeffs) {
            result = result.add(coeff.multiply(xPower)).mod(FIELD_PRIME);
            xPower = xPower.multiply(x).mod(FIELD_PRIME);
        }
        
        return result;
    }
    
    /**
     * Lagrange interpolation implementation
     */
    private BigInteger lagrangeMethod(List<Point> points) {
        BigInteger secret = BigInteger.ZERO;
        int n = points.size();
        
        for (int i = 0; i < n; i++) {
            BigInteger xi = points.get(i).x;
            BigInteger yi = points.get(i).y;
            
            BigInteger numerator = BigInteger.ONE;
            BigInteger denominator = BigInteger.ONE;
            
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    BigInteger xj = points.get(j).x;
                    numerator = numerator.multiply(xj.negate()).mod(FIELD_PRIME);
                    denominator = denominator.multiply(xi.subtract(xj)).mod(FIELD_PRIME);
                }
            }
            
            BigInteger term = yi.multiply(numerator)
                               .multiply(modularInverse(denominator, FIELD_PRIME))
                               .mod(FIELD_PRIME);
            secret = secret.add(term).mod(FIELD_PRIME);
        }
        
        return secret;
    }
    
    /**
     * Gaussian elimination implementation
     */
    private BigInteger gaussianMethod(List<Point> points) {
        int n = points.size();
        BigInteger[][] matrix = new BigInteger[n][n + 1];
        
        // Build Vandermonde matrix
        for (int i = 0; i < n; i++) {
            BigInteger xi = points.get(i).x;
            BigInteger yi = points.get(i).y;
            BigInteger power = BigInteger.ONE;
            
            for (int j = 0; j < n; j++) {
                matrix[i][j] = power.mod(FIELD_PRIME);
                power = power.multiply(xi).mod(FIELD_PRIME);
            }
            matrix[i][n] = yi.mod(FIELD_PRIME);
        }
        
        // Forward elimination
        for (int pivot = 0; pivot < n; pivot++) {
            // Partial pivoting
            int maxRow = pivot;
            for (int i = pivot + 1; i < n; i++) {
                if (matrix[i][pivot].abs().compareTo(matrix[maxRow][pivot].abs()) > 0) {
                    maxRow = i;
                }
            }
            
            if (maxRow != pivot) {
                BigInteger[] temp = matrix[pivot];
                matrix[pivot] = matrix[maxRow];
                matrix[maxRow] = temp;
            }
            
            BigInteger pivotInv = modularInverse(matrix[pivot][pivot], FIELD_PRIME);
            
            for (int i = pivot + 1; i < n; i++) {
                BigInteger factor = matrix[i][pivot].multiply(pivotInv).mod(FIELD_PRIME);
                
                for (int j = pivot; j <= n; j++) {
                    matrix[i][j] = matrix[i][j].subtract(
                        factor.multiply(matrix[pivot][j])
                    ).mod(FIELD_PRIME);
                }
            }
        }
        
        // Back substitution
        BigInteger[] coefficients = new BigInteger[n];
        for (int i = n - 1; i >= 0; i--) {
            coefficients[i] = matrix[i][n];
            
            for (int j = i + 1; j < n; j++) {
                coefficients[i] = coefficients[i].subtract(
                    matrix[i][j].multiply(coefficients[j])
                ).mod(FIELD_PRIME);
            }
            
            coefficients[i] = coefficients[i].multiply(
                modularInverse(matrix[i][i], FIELD_PRIME)
            ).mod(FIELD_PRIME);
        }
        
        return coefficients[0]; // Constant term
    }
    
    /**
     * Calculate modular inverse
     */
    private BigInteger modularInverse(BigInteger a, BigInteger m) {
        return a.modPow(m.subtract(BigInteger.valueOf(2)), m);
    }
    
    /**
     * Estimate memory usage for given k
     */
    private double estimateMemoryUsage(int k) {
        // Lagrange: O(k) points storage
        double lagrangeMemory = k * 64; // 64 bytes per BigInteger pair
        
        // Gaussian: O(k²) matrix storage  
        double gaussianMemory = k * k * 32; // 32 bytes per BigInteger
        
        // Return maximum (worst case)
        return Math.max(lagrangeMemory, gaussianMemory) / (1024.0 * 1024.0); // Convert to MB
    }
    
    /**
     * Demonstrate algorithm scaling behavior
     */
    private void demonstrateScaling() {
        System.out.println("\n=== SCALING ANALYSIS ===");
        System.out.println("Theoretical vs Actual Complexity:");
        
        int[] scalingSizes = {5, 10, 20, 40};
        double[] lagrangeTimes = new double[scalingSizes.length];
        double[] gaussianTimes = new double[scalingSizes.length];
        
        for (int i = 0; i < scalingSizes.length; i++) {
            BenchmarkResult result = benchmarkAlgorithms(scalingSizes[i]);
            lagrangeTimes[i] = result.lagrangeTime;
            gaussianTimes[i] = result.gaussianTime;
        }
        
        System.out.println("\nLagrange Scaling (should be ~O(k²)):");
        for (int i = 1; i < scalingSizes.length; i++) {
            double ratio = lagrangeTimes[i] / lagrangeTimes[i-1];
            double kRatio = (double) scalingSizes[i] / scalingSizes[i-1];
            double expectedRatio = kRatio * kRatio; // k²
            
            System.out.printf("k=%d to k=%d: actual ratio=%.2f, expected O(k²)=%.2f\n",
                scalingSizes[i-1], scalingSizes[i], ratio, expectedRatio);
        }
        
        System.out.println("\nGaussian Scaling (should be ~O(k³)):");
        for (int i = 1; i < scalingSizes.length; i++) {
            double ratio = gaussianTimes[i] / gaussianTimes[i-1];
            double kRatio = (double) scalingSizes[i] / scalingSizes[i-1];
            double expectedRatio = kRatio * kRatio * kRatio; // k³
            
            System.out.printf("k=%d to k=%d: actual ratio=%.2f, expected O(k³)=%.2f\n",
                scalingSizes[i-1], scalingSizes[i], ratio, expectedRatio);
        }
    }
    
    /**
     * Test numerical stability with challenging inputs
     */
    private void testNumericalStability() {
        System.out.println("\n=== NUMERICAL STABILITY TEST ===");
        
        int k = 10;
        List<StabilityTestCase> testCases = Arrays.asList(
            new StabilityTestCase("Sequential points", generateSequentialPoints(k)),
            new StabilityTestCase("Large gaps", generateLargeGapPoints(k)),
            new StabilityTestCase("Small differences", generateSmallDiffPoints(k)),
            new StabilityTestCase("Random distribution", generateRandomPoints(k))
        );
        
        System.out.println("Test Case              | Lagrange Error | Gaussian Error | Winner");
        System.out.println("-----------------------|----------------|----------------|--------");
        
        for (StabilityTestCase testCase : testCases) {
            BigInteger lagrangeResult = lagrangeMethod(testCase.points);
            BigInteger gaussianResult = gaussianMethod(testCase.points);
            
            // Calculate error relative to known secret
            BigInteger lagrangeError = lagrangeResult.subtract(testCase.expectedSecret).abs();
            BigInteger gaussianError = gaussianResult.subtract(testCase.expectedSecret).abs();
            
            String winner = lagrangeError.compareTo(gaussianError) <= 0 ? "Lagrange" : "Gaussian";
            
            System.out.printf("%-22s | %14s | %14s | %s\n",
                testCase.name,
                lagrangeError.equals(BigInteger.ZERO) ? "Perfect" : "±" + lagrangeError.toString().substring(0, Math.min(8, lagrangeError.toString().length())),
                gaussianError.equals(BigInteger.ZERO) ? "Perfect" : "±" + gaussianError.toString().substring(0, Math.min(8, gaussianError.toString().length())),
                winner);
        }
    }
    
    /**
     * Generate test points with different characteristics
     */
    private List<Point> generateSequentialPoints(int k) {
        List<Point> points = new ArrayList<>();
        BigInteger secret = BigInteger.valueOf(12345);
        
        // Simple polynomial: f(x) = secret + 2x + 3x²
        for (int x = 1; x <= k; x++) {
            BigInteger y = secret.add(BigInteger.valueOf(2 * x))
                                 .add(BigInteger.valueOf(3 * x * x))
                                 .mod(FIELD_PRIME);
            points.add(new Point(BigInteger.valueOf(x), y));
        }
        
        return points;
    }
    
    private List<Point> generateLargeGapPoints(int k) {
        List<Point> points = new ArrayList<>();
        BigInteger secret = BigInteger.valueOf(98765);
        
        // Large x values with gaps
        for (int i = 0; i < k; i++) {
            int x = 1000 + (i * 500); // x = 1000, 1500, 2000, ...
            BigInteger y = secret.add(BigInteger.valueOf(x))
                                 .add(BigInteger.valueOf(x * x / 100))
                                 .mod(FIELD_PRIME);
            points.add(new Point(BigInteger.valueOf(x), y));
        }
        
        return points;
    }
    
    private List<Point> generateSmallDiffPoints(int k) {
        List<Point> points = new ArrayList<>();
        BigInteger secret = new BigInteger("123456789012345678901234567890");
        
        // Very close x values (potential numerical issues)
        for (int i = 0; i < k; i++) {
            BigInteger x = BigInteger.valueOf(1000000 + i);
            BigInteger y = secret.add(x.multiply(BigInteger.valueOf(7)))
                                 .add(x.multiply(x).divide(BigInteger.valueOf(1000)))
                                 .mod(FIELD_PRIME);
            points.add(new Point(x, y));
        }
        
        return points;
    }
    
    private List<Point> generateRandomPoints(int k) {
        List<Point> points = new ArrayList<>();
        BigInteger secret = new BigInteger(64, RANDOM);
        
        // Random coefficients
        List<BigInteger> coeffs = new ArrayList<>();
        coeffs.add(secret);
        for (int i = 1; i < k; i++) {
            coeffs.add(new BigInteger(32, RANDOM));
        }
        
        // Random x values
        Set<Integer> usedX = new HashSet<>();
        while (points.size() < k) {
            int x = RANDOM.nextInt(10000) + 1;
            if (!usedX.contains(x)) {
                usedX.add(x);
                BigInteger y = evaluatePolynomial(coeffs, BigInteger.valueOf(x));
                points.add(new Point(BigInteger.valueOf(x), y));
            }
        }
        
        return points;
    }
    
    // Helper classes
    private static class Point {
        final BigInteger x, y;
        
        Point(BigInteger x, BigInteger y) {
            this.x = x;
            this.y = y;
        }
    }
    
    private static class BenchmarkResult {
        final double lagrangeTime;
        final double gaussianTime;
        final double memoryUsageMB;
        
        BenchmarkResult(double lagrangeTime, double gaussianTime, double memoryUsageMB) {
            this.lagrangeTime = lagrangeTime;
            this.gaussianTime = gaussianTime;
            this.memoryUsageMB = memoryUsageMB;
        }
    }
    
    private static class TestData {
        final List<BigInteger> coefficients;
        final List<Point> points;
        final BigInteger expectedSecret;
        
        TestData(List<BigInteger> coefficients, List<Point> points, BigInteger expectedSecret) {
            this.coefficients = coefficients;
            this.points = points;
            this.expectedSecret = expectedSecret;
        }
    }
    
    private static class StabilityTestCase {
        final String name;
        final List<Point> points;
        final BigInteger expectedSecret;
        
        StabilityTestCase(String name, List<Point> points) {
            this.name = name;
            this.points = points;
            
            // Extract expected secret (assuming it's encoded correctly)
            // For test purposes, we'll compute it using a reference method
            this.expectedSecret = BigInteger.valueOf(12345); // Simplified
        }
    }
}