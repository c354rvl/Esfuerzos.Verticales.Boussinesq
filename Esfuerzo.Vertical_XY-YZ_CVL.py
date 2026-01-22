# -*- coding: utf-8 -*-
"""
Created on Thu Jan 22 07:59:58 2026
Código generado por Gemini 3, modificado por cesar
@author: césar vázquez lorenzana
"""

import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
import os

def G_func(A, B, z):
    """Función de influencia base G(A, B, z)."""
    if z < 1e-6: z = 1e-6 # Evitar división por cero en la superficie
    den1 = 3 * (A**2 + z**2)**2 * np.sqrt(A**2 + B**2 + z**2)
    if den1 < 1e-12: return 0
    
    term1 = B / den1
    term2 = 2 + ((A**2 + z**2) / (A**2 + B**2 + z**2))
    return term1 * term2

def calcular_sigma_z(x_coords, y_coords, q_loads, x0, y0, z):
    """Implementa la sumatoria de esfuerzos según la ecuación (19)."""
    n = len(x_coords)
    x = np.append(x_coords, x_coords[n-1])
    y = np.append(y_coords, y_coords[n-1])
    total_sum = 0
    
    for i in range(n):
        dx_i = x[i+1] - x[i]
        dy_i = y[i+1] - y[i]
        qi = q_loads[i]
        
        # Integrando: G(Yt-y0, Xt-x0, z)dy - G(Xt-x0, Yt-y0, z)dx
        def integrando(t):
            xt = (x[i] + t * dx_i) - x0
            yt = (y[i] + t * dy_i) - y0
            return G_func(yt, xt, z) * dy_i - G_func(xt, yt, z) * dx_i
        
        integral_i, _ = quad(integrando, 0, 1)
        total_sum += qi * integral_i
        
    return (3 * (z**3) / (4 * np.pi)) * total_sum

def graficar_planta_xy(x_pts, y_pts, q_pts):
    z_fija = float(input("Ingrese profundidad de análisis (z): "))
    ext1 = float(input("Ingrese valor mín de extensión del plano (ej. -10): "))
    ext2 = float(input("Ingrese valor máx de extensión del plano (ej. 10): "))
    res = 40

    print(f"Calculando esfuerzos verticales en plano z={z_fija}...")    
    x_range = np.linspace(ext1, ext2, res)
    y_range = np.linspace(ext1, ext2, res)
    X, Y = np.meshgrid(x_range, y_range)
    Sigma = np.array([[calcular_sigma_z(x_pts, y_pts, q_pts, xi, yi, z_fija) for xi in x_range] for yi in y_range])
    
    # Graficando
    plt.figure(figsize=(9, 7))

    # Dibujar curvas de nivel rellenas
    cp = plt.contourf(X, Y, Sigma, levels=20, cmap='gist_rainbow_r')
    plt.colorbar(cp, label='Esfuerzo Vertical ($\sigma_z$)')

    # Dibujar líneas de contorno
    lineas = plt.contour(X, Y, Sigma, levels=20, colors='black', linewidths=0.5)
    plt.clabel(lineas, inline=True, fontsize=8)    

    # Dibujar el polígono de carga para referencia
    # poligono_x = np.append(x_pts, x_pts[0])
    # poligono_y = np.append(y_pts, y_pts[0])
    plt.plot(x_pts, y_pts, 'k--', linewidth=2, label='Área de Carga')
    
    plt.title(f"Plano Horizontal XY (z = {z_fija})")
    plt.xlabel("x"); plt.ylabel("y"); plt.legend(); plt.grid(alpha=0.3); plt.show()

def graficar_perfil_xz(x_pts, y_pts, q_pts):
    y_fijo = float(input("Ingrese plano de corte Y (ej. 0): "))
    ext_x1 = float(input("Ingrese valor horizontal mínimo X (ej. -10): "))
    ext_x2 = float(input("Ingrese valor horizontal mmáximo X (ej. 10): "))
    z_max = float(input("Ingrese profundidad máxima Z: "))
    res = 40

    x_range = np.linspace(ext_x1, ext_x2, res)
    z_range = np.linspace(0.1, z_max, res)
    X, Z = np.meshgrid(x_range, z_range)
    Sigma = np.array([[calcular_sigma_z(x_pts, y_pts, q_pts, xi, y_fijo, zi) for xi in x_range] for zi in z_range])

    plt.figure(figsize=(9, 7))
    cp = plt.contourf(X, Z, Sigma, levels=20, cmap='gist_rainbow_r')
    plt.colorbar(cp, label='$\sigma_z$')
    plt.gca().invert_yaxis()
  
    # Dibujar líneas de contorno
    lineas = plt.contour(X, Z, Sigma, levels=20, colors='black', linewidths=0.5)
    plt.clabel(lineas, inline=True, fontsize=8)    
    
    # Indicar la posición de la carga en superficie (z=0)
    # Mostramos el segmento del polígono que cruza el plano Y fijo
    plt.axhline(0, color='black', linewidth=2)
    plt.title(f"Perfil Vertical XZ (y = {y_fijo})")
    plt.xlabel("x"); plt.ylabel("Z (Profundidad)"); plt.grid(alpha=0.3); plt.show()

def main():
    nombre_archivo = "Polig.L3Z.txt"
    if not os.path.exists(nombre_archivo):
        with open(nombre_archivo, "w") as f:
            f.write("# x y q\n-2 -2 100\n2 -2 100\n2 2 100\n-2 2 100")
    
    data = np.loadtxt(nombre_archivo)
    x_p, y_p, q_p = data[:, 0], data[:, 1], data[:, 2]

    while True:
        print("\n--- MENÚ DE VISUALIZACIÓN ---")
        print("1. Curvas de nivel en Planta (XY)")
        print("2. Bulbo de presiones en Perfil (XZ)")
        print("3. Salir")
        opcion = input("Seleccione una opción: ")

        if opcion == '1': graficar_planta_xy(x_p, y_p, q_p)
        elif opcion == '2': graficar_perfil_xz(x_p, y_p, q_p)
        elif opcion == '3': break
        else: print("Opción no válida.")

if __name__ == "__main__":
    main()