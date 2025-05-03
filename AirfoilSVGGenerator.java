import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class AirfoilSVGGenerator {

    // Cubic spline with support for context-aware interpolation
    public static class Spline {
	private final double[] x, y, a, b, c, d;
    
	public Spline(double[] x, double[] y) {
	    int n = x.length;
	    this.x = x;
	    this.y = y;
	    a = new double[n];
	    b = new double[n];
	    c = new double[n];
	    d = new double[n];
	    
	    double[] h = new double[n - 1];
	    double[] alpha = new double[n - 1];
	    for (int i = 0; i < n - 1; i++) {
		h[i] = x[i + 1] - x[i];
		alpha[i] = (y[i + 1] - y[i]) / h[i];
	    }
	    
	    double[] l = new double[n];
	    double[] mu = new double[n];
	    double[] z = new double[n];
	
	    l[0] = 1;
	    for (int i = 1; i < n - 1; i++) {
		l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
		mu[i] = h[i] / l[i];
		z[i] = (3 * (y[i + 1] - y[i]) / h[i] - 3 * (y[i] - y[i - 1]) / h[i - 1] - mu[i - 1] * z[i - 1]) / l[i];
	    }
	    l[n - 1] = 1;
	    
	    for (int j = n - 2; j >= 0; j--) {
		c[j] = z[j] - mu[j] * c[j + 1];
		b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
		d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
		a[j] = y[j];
	    }
	}

	public double interpolate(double xi) {
	    int i = Arrays.binarySearch(x, xi);
	    if (i < 0) i = -i - 2;
	    if (i >= x.length - 1) i = x.length - 2;
	    double dx = xi - x[i];
	    return a[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
	}
    }

    // Build a smooth curve for each 5% segment using nearby points
    public static Spline buildContextAwareSpline(double[] x, double[] y, double start, double end, double range) {
	List<Double> ctxX = new ArrayList<>();
	List<Double> ctxY = new ArrayList<>();
	for (int i = 0; i < x.length; i++) {
	    if (x[i] >= (start - range) && x[i] <= (end + range)) {
		ctxX.add(x[i]);
		ctxY.add(y[i]);
	    }
	}
	double[] cx = ctxX.stream().mapToDouble(Double::doubleValue).toArray();
	double[] cy = ctxY.stream().mapToDouble(Double::doubleValue).toArray();
	return new Spline(cx, cy);
    }
    
    public static void main(String[] args) throws IOException {
	double[] chordPercent = new double[]{0.0, 0.025, 0.05, 0.075, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 1.0};
	double[] upperY = new double[]{0.035, 0.065, 0.079, 0.0885, 0.096, 0.1069, 0.1136, 0.1170, 0.1140, 0.1052, 0.0915, 0.0735, 0.0522, 0.0280, 0.0149, 0.0012};
	double[] lowerY = new double[]{0.035,0.0147,0.0093,0.0063,0.0042,0.0015,0.0003, 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
	
    StringBuilder path = new StringBuilder("M ");

    for (double x = 0.0; x <= 1.0; x += 0.01) {
        Spline upperSpline = buildContextAwareSpline(chordPercent, upperY, x - 0.05, x + 0.05, 0.05);
        double y = upperSpline.interpolate(x);
        path.append(String.format("%.4f,%.4f ", x * 1000, (1 - y) * 1000));
    }
    for (double x = 1.0; x >= 0.0; x -= 0.01) {
        Spline lowerSpline = buildContextAwareSpline(chordPercent, lowerY, x - 0.05, x + 0.05, 0.05);
        double y = lowerSpline.interpolate(x);
        path.append(String.format("%.4f,%.4f ", x * 1000, (1 - y) * 1000));
    }
    path.append("Z");

    try (FileWriter writer = new FileWriter("clark_y_context.svg")) {
        writer.write("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1000\" height=\"1000\" viewBox=\"0 0 1000 1000\">\n");
        writer.write("<path d=\"" + path.toString() + "\" fill=\"none\" stroke=\"black\"/>\n");
        writer.write("</svg>");
    }

    System.out.println("SVG with context-aware spline saved to clark_y_context.svg");
    }

}
