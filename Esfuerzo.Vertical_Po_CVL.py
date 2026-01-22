# -*- coding: utf-8 -*-
"""
Created on Thu Jan 22 07:25:52 2026
Código generado por Gemini 3, modificado por cesar
@author: césar vázquez lorenzana
"""
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
import os

def G_func(A, B, z):
    """Función de influencia base G(A, B, z)."""
    den1 = 3 * (A**2 + z**2)**2 * np.sqrt(A**2 + B**2 + z**2)
    if den1 < 1e-12: # Evitar división por cero en bordes
        return 0
    
    term1 = B / den1
    term2 = 2 + ((A**2 + z**2) / (A**2 + B**2 + z**2))
    return term1 * term2

def calcular_sigma_z(x_coords, y_coords, q_loads, x0, y0, z):
    """
    Calcula el esfuerzo vertical sigma_z en el punto P(x0, y0, z).
    Implementa la ecuación (19) de la imagen proporcionada.
    """
    n = len(x_coords)
    # Cerrar el polígono conectando el último punto con el primero
    x = np.append(x_coords, x_coords[n-1])
    y = np.append(y_coords, y_coords[n-1])
    
    total_sum = 0
    
    for i in range(n):
        dx_i = x[i+1] - x[i] # Corresponde a Delta x_i
        dy_i = y[i+1] - y[i] # Corresponde a Delta y_i
        qi = q_loads[i]      # Carga del segmento i
        
        # Funciones paramétricas trasladadas al punto (x0, y0)
        def Xi_trasladada(t): return (x[i] + t * dx_i) - x0
        def Yi_trasladada(t): return (y[i] + t * dy_i) - y0
        
        # Integrando de la ecuación (19)
        def integrando(t):
            term_a = G_func(Yi_trasladada(t), Xi_trasladada(t), z) * dy_i
            term_b = G_func(Xi_trasladada(t), Yi_trasladada(t), z) * dx_i
            return term_a - term_b
        
        integral_i, _ = quad(integrando, 0, 1)
        total_sum += qi * integral_i
        
    return (3 * (z**3) / (4 * np.pi)) * total_sum

def main():
    nombre_archivo = "Polig.L3R.txt"
    
    # Crear archivo de ejemplo si no existe (Formato: x y q)
    if not os.path.exists(nombre_archivo):
        with open(nombre_archivo, "w") as f:
            f.write("# x y q\n-1 -1 100\n1 -1 100\n1 1 100\n-1 1 100")
        print(f"Archivo '{nombre_archivo}' creado con datos de ejemplo.")

    try:
        # Carga x, y, q separadlos por espacios
        data = np.loadtxt(nombre_archivo)
        x_pts, y_pts, q_pts = data[:, 0], data[:, 1], data[:, 2]
        
        print("\n--- Configuración del Punto de Observación ---")
        x0 = float(input("Ingrese coordenada x0: "))
        y0 = float(input("Ingrese coordenada y0: "))
        z  = float(input("Ingrese profundidad z: "))
        
        resultado = calcular_sigma_z(x_pts, y_pts, q_pts, x0, y0, z)

        print(f"Archivo '{nombre_archivo}'")
        # print(f"\nResultado:")
        print(f"Esfuerzo vertical (sigma_z) en ({x0}, {y0}, {z}) = {resultado:.6f}")
        # Graficando
        plt.figure(figsize=(9, 7))
        # Dibujar el polígono de carga para referencia
        plt.plot(x_pts, y_pts, 'k--', linewidth=2, label='Área de Carga')
        plt.title(f"Plano Horizontal XY ")
        plt.xlabel("x"); plt.ylabel("y"); plt.legend(); plt.grid(alpha=0.3); plt.show()
        
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()