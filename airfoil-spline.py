
import numpy as np

# Spline interpolator (natural cubic) without external libraries
def natural_cubic_spline(x, y):
    n = len(x)
    a = y[:]
    h = [x[i + 1] - x[i] for i in range(n - 1)]
    alpha = [0] * n
    for i in range(1, n - 1):
        alpha[i] = (3 / h[i]) * (a[i + 1] - a[i]) - (3 / h[i - 1]) * (a[i] - a[i - 1])

    l = [1] + [0] * (n - 1)
    mu = [0] * n
    z = [0] * n
    for i in range(1, n - 1):
        l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1]
        mu[i] = h[i] / l[i]
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i]
    l[-1] = 1

    c = [0] * n
    b = [0] * (n - 1)
    d = [0] * (n - 1)
    for j in range(n - 2, -1, -1):
        c[j] = z[j] - mu[j] * c[j + 1]
        b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3
        d[j] = (c[j + 1] - c[j]) / (3 * h[j])
    return a, b, c, d, x

def evaluate_spline(a, b, c, d, x_orig, x_eval):
    result = []
    for x in x_eval:
        i = np.searchsorted(x_orig, x) - 1
        if i < 0:
            i = 0
        elif i >= len(b):
            i = len(b) - 1
        dx = x - x_orig[i]
        y = a[i] + b[i] * dx + c[i] * dx ** 2 + d[i] * dx ** 3
        result.append(y)
    return result

# Data
chord = [0.0, 0.025, 0.05, 0.075, 0.10, 0.15, 0.20,
         0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 1.0]
upper = [0.035, 0.065, 0.079, 0.0885, 0.096, 0.1069, 0.1136,
         0.1170, 0.1140, 0.1052, 0.0915, 0.0735, 0.0522, 0.0280, 0.0149, 0.0012]
lower = [0.035, 0.0147, 0.0093, 0.0063, 0.0042, 0.0015, 0.0003,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# Spline calculation
a_u, b_u, c_u, d_u, x_u = natural_cubic_spline(chord, upper)
a_l, b_l, c_l, d_l, x_l = natural_cubic_spline(chord, lower)

x_dense = np.linspace(0, 1, 200)
y_u_dense = evaluate_spline(a_u, b_u, c_u, d_u, x_u, x_dense)
y_l_dense = evaluate_spline(a_l, b_l, c_l, d_l, x_l, x_dense)

# SVG path generation
path_data = 'M '
for x, y in zip(x_dense, y_u_dense):
    path_data += f'{x*1000:.2f},{(1 - y)*1000:.2f} '
for x, y in zip(reversed(x_dense), reversed(y_l_dense)):
    path_data += f'{x*1000:.2f},{(1 - y)*1000:.2f} '
path_data += 'Z'

# Save SVG file to the same directory
svg_content = f'''<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="1000" viewBox="0 0 1000 1000">
<path d="{path_data}" fill="none" stroke="black"/>
</svg>
'''

with open("airfoil_spline22.svg", "w") as f:
    f.write(svg_content)

print("SVG file 'airfoil_spline.svg' has been created in the current directory.")
