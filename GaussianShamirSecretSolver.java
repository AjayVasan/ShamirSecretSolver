import java.io.*;
import java.math.BigInteger;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.IntStream;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;

/**
 * Enhanced Shamir's Secret Sharing Solver
 * Uses Gaussian Elimination for large k, Lagrange for small k
 * Production-grade implementation with modular arithmetic optimization
 */
public class GaussianShamirSecretSolver {
    
    // Large prime for finite field arithmetic (2^127 - 1, Mersenne prime)
    private static final BigInteger FIELD_PRIME = new BigInteger("170141183460469231731687303715884105727");
    
    // Threshold for algorithm selection
    private static final int GAUSSIAN_THRESHOLD = 10;
    
    // Thread pool for parallel processing
    private static final ForkJoinPool THREAD_POOL = new ForkJoinPool();
    
    public static void main(String[] args) {
        ShamirSecretSolver solver = new ShamirSecretSolver();
        
        try {
            System.out.println("=== ENHANCED SHAMIR'S SECRET SHARING SOLVER ===");
            System.out.println("Algorithm Selection: Lagrange (k≤10) | Gaussian Elimination (k>10)");
            System.out.println("Field Prime: " + FIELD_PRIME.toString().substring(0, 20) + "...\n");
            
            // Process both test cases with performance metrics
            long startTime = System.nanoTime();
            
            // Test Case 1
            System.out.println("Processing Test Case 1:");
            BigInteger secret1 = solver.solveFromFile("testcase1.json");
            System.out.println("Secret for Test Case 1: " + secret1);
            System.out.println();
            
            // Test Case 2
            System.out.println("Processing Test Case 2:");
            BigInteger secret2 = solver.solveFromFile("testcase2.json");
            System.out.println("Secret for Test Case 2: " + secret2);
            System.out.println();
            
            long endTime = System.nanoTime();
            double executionTime = (endTime - startTime) / 1_000_000.0; // Convert to milliseconds
            
            // Final Results
            System.out.println("=== FINAL RESULTS ===");
            System.out.println("Test Case 1 Secret: " + secret1);
            System.out.println("Test Case 2 Secret: " + secret2);
            System.out.println("Total Execution Time: " + String.format("%.2f ms", executionTime));
            
        } catch (Exception e) {
            System.err.println("Error: " + e.getMessage());
            e.printStackTrace();
        } finally {
            THREAD_POOL.shutdown();
        }
    }
    
    /**
     * Solve secret from JSON file
     */
    public BigInteger solveFromFile(String filename) throws IOException {
        ObjectMapper mapper = new ObjectMapper();
        JsonNode root = mapper.readTree(new File(filename));
        return solveSecret(root);
    }
    
    /**
     * Enhanced solver with algorithm selection
     */
    public BigInteger solveSecret(JsonNode root) {
        // Parse input parameters
        JsonNode keysNode = root.get("keys");
        int n = keysNode.get("n").asInt();
        int k = keysNode.get("k").asInt();
        
        System.out.println("n (total roots): " + n);
        System.out.println("k (required roots): " + k);
        System.out.println("Polynomial degree: " + (k - 1));
        
        // Parse and decode all roots
        List<Point> allPoints = parseAndDecodeRoots(root);
        
        // Select algorithm based on problem size
        String algorithm = (k <= GAUSSIAN_THRESHOLD) ? "Lagrange Interpolation" : "Gaussian Elimination";
        System.out.println("Selected Algorithm: " + algorithm);
        System.out.println("Complexity: O(" + (k <= GAUSSIAN_THRESHOLD ? "k²" : "k³") + ")");
        
        // Select first k points for computation
        List<Point> selectedPoints = allPoints.subList(0, Math.min(k, allPoints.size()));
        
        System.out.println("\nUsing " + selectedPoints.size() + " points for computation:");
        for (int i = 0; i < Math.min(5, selectedPoints.size()); i++) {
            Point p = selectedPoints.get(i);
            System.out.println("(" + p.x + ", " + p.y.toString().substring(0, Math.min(15, p.y.toString().length())) + "...)");
        }
        if (selectedPoints.size() > 5) {
            System.out.println("... and " + (selectedPoints.size() - 5) + " more points");
        }
        
        // Apply selected algorithm
        BigInteger secret;
        long algorithmStart = System.nanoTime();
        
        if (k <= GAUSSIAN_THRESHOLD) {
            secret = lagrangeInterpolation(selectedPoints);
        } else {
            secret = gaussianEliminationMethod(selectedPoints);
        }
        
        long algorithmEnd = System.nanoTime();
        double algorithmTime = (algorithmEnd - algorithmStart) / 1_000_000.0;
        
        System.out.println("Algorithm execution time: " + String.format("%.3f ms", algorithmTime));
        
        // Verification with cross-validation
        if (allPoints.size() > k) {
            verifyWithCrossValidation(allPoints, k, secret);
        }
        
        return secret;
    }
    
    /**
     * Parse and decode roots from JSON with enhanced error handling
     */
    private List<Point> parseAndDecodeRoots(JsonNode root) {
        List<Point> allPoints = new ArrayList<>();
        
        System.out.println("\nDecoding roots:");
        Iterator<String> fieldNames = root.fieldNames();
        
        while (fieldNames.hasNext()) {
            String fieldName = fieldNames.next();
            
            if ("keys".equals(fieldName)) continue;
            
            try {
                int x = Integer.parseInt(fieldName);
                JsonNode pointNode = root.get(fieldName);
                
                int base = pointNode.get("base").asInt();
                String valueStr = pointNode.get("value").asText();
                
                // Enhanced base validation
                if (base < 2 || base > 36) {
                    throw new IllegalArgumentException("Invalid base: " + base + " (must be 2-36)");
                }
                
                // Decode with overflow protection
                BigInteger y = decodeValueSafely(valueStr, base);
                allPoints.add(new Point(BigInteger.valueOf(x), y));
                
                String displayValue = valueStr.length() > 20 ? 
                    valueStr.substring(0, 20) + "..." : valueStr;
                String displayDecoded = y.toString().length() > 20 ? 
                    y.toString().substring(0, 20) + "..." : y.toString();
                
                System.out.println("Root " + x + ": base=" + base + 
                                 ", encoded='" + displayValue + 
                                 "', decoded=" + displayDecoded);
                
            } catch (NumberFormatException e) {
                continue; // Skip non-numeric field names
            } catch (Exception e) {
                System.err.println("Warning: Error processing field " + fieldName + ": " + e.getMessage());
            }
        }
        
        // Sort points by x coordinate for better numerical stability
        allPoints.sort(Comparator.comparing(p -> p.x));
        
        return allPoints;
    }
    
    /**
     * Safe base conversion with overflow detection
     */
    private BigInteger decodeValueSafely(String value, int base) {
        try {
            // Pre-validate characters for the given base
            for (char c : value.toCharArray()) {
                int digit = Character.digit(c, base);
                if (digit == -1) {
                    throw new IllegalArgumentException("Invalid character '" + c + "' for base " + base);
                }
            }
            
            BigInteger result = new BigInteger(value, base);
            
            // Ensure result fits in our field
            if (result.compareTo(FIELD_PRIME) >= 0) {
                result = result.mod(FIELD_PRIME);
            }
            
            return result;
            
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Invalid value '" + value + "' for base " + base, e);
        }
    }
    
    /**
     * Enhanced Lagrange interpolation with modular arithmetic
     */
    private BigInteger lagrangeInterpolation(List<Point> points) {
        System.out.println("\nApplying Enhanced Lagrange Interpolation:");
        
        BigInteger secret = BigInteger.ZERO;
        int n = points.size();
        
        for (int i = 0; i < n; i++) {
            BigInteger xi = points.get(i).x;
            BigInteger yi = points.get(i).y;
            
            // Calculate Lagrange basis polynomial Li(0) with lazy reduction
            BigInteger numerator = BigInteger.ONE;
            BigInteger denominator = BigInteger.ONE;
            
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    BigInteger xj = points.get(j).x;
                    
                    // numerator *= (0 - xj) = -xj
                    numerator = numerator.multiply(xj.negate());
                    
                    // denominator *= (xi - xj)
                    denominator = denominator.multiply(xi.subtract(xj));
                }
            }
            
            // Modular arithmetic to avoid division
            BigInteger denominatorInv = modularInverse(denominator, FIELD_PRIME);
            BigInteger term = yi.multiply(numerator).multiply(denominatorInv).mod(FIELD_PRIME);
            
            secret = secret.add(term).mod(FIELD_PRIME);
            
            if (i < 3 || n <= 5) { // Show details for small cases
                System.out.println("  L" + i + "(0) contribution: " + 
                                 yi.toString().substring(0, Math.min(10, yi.toString().length())) + "... * (" + 
                                 numerator.toString().substring(0, Math.min(8, numerator.toString().length())) + "... / " + 
                                 denominator.toString().substring(0, Math.min(8, denominator.toString().length())) + "...) = " + 
                                 term.toString().substring(0, Math.min(10, term.toString().length())) + "...");
            }
        }
        
        System.out.println("Lagrange secret: " + secret);
        return secret;
    }
    
    /**
     * Gaussian Elimination method with finite field arithmetic
     */
    private BigInteger gaussianEliminationMethod(List<Point> points) {
        System.out.println("\nApplying Gaussian Elimination over GF(p):");
        
        int n = points.size();
        
        // Build Vandermonde matrix [x^0, x^1, x^2, ..., x^(n-1)]
        // System: Ax = y, where A[i][j] = x_i^j
        BigInteger[][] matrix = new BigInteger[n][n + 1]; // Augmented matrix
        
        System.out.println("Building Vandermonde matrix...");
        
        // Parallel matrix construction for large k
        if (n > 20) {
            THREAD_POOL.submit(() -> {
                IntStream.range(0, n).parallel().forEach(i -> {
                    BigInteger xi = points.get(i).x;
                    BigInteger yi = points.get(i).y;
                    BigInteger power = BigInteger.ONE;
                    
                    for (int j = 0; j < n; j++) {
                        matrix[i][j] = power.mod(FIELD_PRIME);
                        power = power.multiply(xi).mod(FIELD_PRIME);
                    }
                    matrix[i][n] = yi.mod(FIELD_PRIME); // Augmented column
                });
            }).join();
        } else {
            // Sequential construction for small k
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
        }
        
        // Forward elimination
        System.out.println("Forward elimination phase...");
        for (int pivot = 0; pivot < n; pivot++) {
            // Find pivot row (partial pivoting for numerical stability)
            int maxRow = pivot;
            for (int i = pivot + 1; i < n; i++) {
                if (matrix[i][pivot].abs().compareTo(matrix[maxRow][pivot].abs()) > 0) {
                    maxRow = i;
                }
            }
            
            // Swap rows if needed
            if (maxRow != pivot) {
                BigInteger[] temp = matrix[pivot];
                matrix[pivot] = matrix[maxRow];
                matrix[maxRow] = temp;
            }
            
            // Check for zero pivot
            if (matrix[pivot][pivot].equals(BigInteger.ZERO)) {
                throw new ArithmeticException("Matrix is singular - cannot solve system");
            }
            
            // Eliminate column entries below pivot
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
        System.out.println("Back substitution phase...");
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
        
        // The constant term is coefficients[0]
        BigInteger secret = coefficients[0];
        
        System.out.println("Polynomial coefficients found:");
        for (int i = 0; i < Math.min(5, n); i++) {
            System.out.println("  a" + i + " = " + 
                             coefficients[i].toString().substring(0, Math.min(15, coefficients[i].toString().length())) + "...");
        }
        if (n > 5) {
            System.out.println("  ... and " + (n - 5) + " more coefficients");
        }
        
        System.out.println("Gaussian elimination secret: " + secret);
        return secret;
    }
    
    /**
     * Extended Euclidean algorithm for modular inverse
     */
    private BigInteger modularInverse(BigInteger a, BigInteger m) {
        if (a.gcd(m).equals(BigInteger.ONE)) {
            return a.modPow(m.subtract(BigInteger.valueOf(2)), m);
        }
        throw new ArithmeticException("Modular inverse does not exist for " + a + " mod " + m);
    }
    
    /**
     * Cross-validation with multiple point combinations
     */
    private void verifyWithCrossValidation(List<Point> allPoints, int k, BigInteger expectedSecret) {
        System.out.println("\nCross-validation:");
        
        int validations = Math.min(3, allPoints.size() - k + 1);
        int successCount = 0;
        
        for (int offset = 1; offset <= validations; offset++) {
            try {
                List<Point> testPoints = allPoints.subList(offset, offset + k);
                BigInteger testSecret;
                
                if (k <= GAUSSIAN_THRESHOLD) {
                    testSecret = lagrangeInterpolation(testPoints);
                } else {
                    testSecret = gaussianEliminationMethod(testPoints);
                }
                
                boolean matches = testSecret.equals(expectedSecret);
                if (matches) successCount++;
                
                System.out.println("  Validation " + offset + ": " + 
                                 (matches ? "✓ PASS" : "✗ FAIL") + 
                                 " (secret: " + testSecret.toString().substring(0, Math.min(15, testSecret.toString().length())) + "...)");
                
            } catch (Exception e) {
                System.out.println("  Validation " + offset + ": ✗ ERROR - " + e.getMessage());
            }
        }
        
        System.out.println("Cross-validation result: " + successCount + "/" + validations + " passed");
        
        if (successCount == validations) {
            System.out.println("✅ All validations passed - high confidence in result");
        } else if (successCount > validations / 2) {
            System.out.println("⚠️ Majority passed - moderate confidence");
        } else {
            System.out.println("❌ Low confidence - possible data corruption");
        }
    }
    
    /**
     * Point class with enhanced functionality
     */
    private static class Point {
        final BigInteger x, y;
        
        Point(BigInteger x, BigInteger y) {
            this.x = x;
            this.y = y;
        }
        
        @Override
        public String toString() {
            return "(" + x + ", " + y.toString().substring(0, Math.min(10, y.toString().length())) + "...)";
        }
        
        @Override
        public boolean equals(Object obj) {
            if (this == obj) return true;
            if (!(obj instanceof Point)) return false;
            Point other = (Point) obj;
            return x.equals(other.x) && y.equals(other.y);
        }
        
        @Override
        public int hashCode() {
            return Objects.hash(x, y);
        }
    }
}