import java.io.*;
import java.math.BigInteger;
import java.util.*;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;

/**
 * Shamir's Secret Sharing Solver
 * Solves for the constant term of a polynomial using given roots
 * encoded in different bases
 */
public class ShamirSecretSolver {
    
    public static void main(String[] args) {
        ShamirSecretSolver solver = new ShamirSecretSolver();
        
        try {
            // Process both test cases
            System.out.println("=== SHAMIR'S SECRET SHARING SOLVER ===\n");
            
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
            
            // Summary
            System.out.println("=== FINAL RESULTS ===");
            System.out.println("Test Case 1 Secret: " + secret1);
            System.out.println("Test Case 2 Secret: " + secret2);
            
        } catch (Exception e) {
            System.err.println("Error: " + e.getMessage());
            e.printStackTrace();
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
     * Main solver method
     */
    public BigInteger solveSecret(JsonNode root) {
        // Step 1: Parse input parameters
        JsonNode keysNode = root.get("keys");
        int n = keysNode.get("n").asInt(); // Total number of roots provided
        int k = keysNode.get("k").asInt(); // Minimum roots required (degree + 1)
        
        System.out.println("n (total roots): " + n);
        System.out.println("k (required roots): " + k);
        System.out.println("Polynomial degree: " + (k - 1));
        
        // Step 2: Parse and decode all roots
        List<Point> allPoints = new ArrayList<>();
        
        System.out.println("\nDecoding roots:");
        Iterator<String> fieldNames = root.fieldNames();
        while (fieldNames.hasNext()) {
            String fieldName = fieldNames.next();
            
            // Skip the "keys" field
            if ("keys".equals(fieldName)) {
                continue;
            }
            
            try {
                int x = Integer.parseInt(fieldName);
                JsonNode pointNode = root.get(fieldName);
                
                int base = pointNode.get("base").asInt();
                String valueStr = pointNode.get("value").asText();
                
                // Decode the y value from the given base
                BigInteger y = decodeValue(valueStr, base);
                
                allPoints.add(new Point(BigInteger.valueOf(x), y));
                
                System.out.println("Root " + x + ": base=" + base + ", encoded='" + 
                                 valueStr + "', decoded=" + y);
                
            } catch (NumberFormatException e) {
                // Skip non-numeric field names
                continue;
            }
        }
        
        // Step 3: Select first k points for interpolation
        // (We could use any k points, but first k is sufficient)
        List<Point> selectedPoints = allPoints.subList(0, Math.min(k, allPoints.size()));
        
        System.out.println("\nUsing " + selectedPoints.size() + " points for interpolation:");
        for (Point p : selectedPoints) {
            System.out.println("(" + p.x + ", " + p.y + ")");
        }
        
        // Step 4: Apply Lagrange interpolation to find secret (constant term)
        BigInteger secret = lagrangeInterpolation(selectedPoints);
        
        // Step 5: Verification (optional) - test with different combinations
        if (allPoints.size() >= k && k > 1) {
            System.out.println("\nVerification with different point combinations:");
            
            // Try with points 2 to k+1 if available
            if (allPoints.size() > k) {
                List<Point> verifyPoints = allPoints.subList(1, k + 1);
                BigInteger verifySecret = lagrangeInterpolation(verifyPoints);
                
                System.out.println("Verification result: " + verifySecret);
                System.out.println("Secrets match: " + secret.equals(verifySecret));
            }
        }
        
        return secret;
    }
    
    /**
     * Decode value from given base to decimal
     */
    private BigInteger decodeValue(String value, int base) {
        if (base < 2 || base > 36) {
            throw new IllegalArgumentException("Base must be between 2 and 36");
        }
        
        try {
            return new BigInteger(value, base);
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Invalid value '" + value + "' for base " + base);
        }
    }
    
    /**
     * Lagrange interpolation to find the constant term (secret)
     * Evaluates the polynomial at x = 0
     */
    private BigInteger lagrangeInterpolation(List<Point> points) {
        BigInteger secret = BigInteger.ZERO;
        int n = points.size();
        
        System.out.println("\nApplying Lagrange Interpolation:");
        
        for (int i = 0; i < n; i++) {
            BigInteger xi = points.get(i).x;
            BigInteger yi = points.get(i).y;
            
            // Calculate Lagrange basis polynomial Li(0)
            // Li(0) = ∏(0 - xj) / ∏(xi - xj) for j ≠ i
            BigInteger numerator = BigInteger.ONE;
            BigInteger denominator = BigInteger.ONE;
            
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    BigInteger xj = points.get(j).x;
                    
                    // For x = 0: numerator *= (0 - xj) = -xj
                    numerator = numerator.multiply(xj.negate());
                    
                    // denominator *= (xi - xj)
                    denominator = denominator.multiply(xi.subtract(xj));
                }
            }
            
            // Calculate Li(0) = numerator / denominator
            // Since we're working with integers, we need to be careful with division
            if (denominator.equals(BigInteger.ZERO)) {
                throw new ArithmeticException("Division by zero in Lagrange interpolation");
            }
            
            // For exact division, we calculate: yi * numerator / denominator
            BigInteger term = yi.multiply(numerator).divide(denominator);
            secret = secret.add(term);
            
            System.out.println("  L" + i + "(0) contribution: " + yi + " * (" + 
                             numerator + " / " + denominator + ") = " + term);
        }
        
        System.out.println("Final secret (constant term): " + secret);
        return secret;
    }
    
    /**
     * Alternative method using matrix approach (for verification)
     * This method solves the system of linear equations to find coefficients
     */
    public BigInteger solveUsingMatrix(List<Point> points) {
        int n = points.size();
        
        // Create Vandermonde matrix and solve for coefficients
        // This is more complex to implement but provides verification
        // For now, we'll stick with Lagrange interpolation as it's more direct
        
        System.out.println("Matrix method not implemented in this version");
        return BigInteger.ZERO;
    }
    
    /**
     * Point class to represent (x, y) coordinates
     */
    private static class Point {
        BigInteger x, y;
        
        Point(BigInteger x, BigInteger y) {
            this.x = x;
            this.y = y;
        }
        
        @Override
        public String toString() {
            return "(" + x + ", " + y + ")";
        }
    }
}