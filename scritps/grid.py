from glob import glob
import subprocess
import re
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import pandas as pd

USER_MENU =\
"""
Category:
 1) Homogeneous = 1
 2) Homogeneous > 1
 3) Cold spots
 4) Hot spots
 5) Stripes
(6) Competition
(7) Chaos
"""

def fill_grid():
    with open('classification.csv', 'w') as f:
        f.write('phi,gamma,category\n')


    for filename in tqdm(glob(r'.\output\frame*.svg')):
        
        print(filename)
        match = re.search(r'frame_u_phi([0-9.]+)_gamma([0-9.]+)_', filename)
        if match:
            phi = match.group(1)
            gamma = match.group(2)
        
        data = np.loadtxt(filename.replace('.svg', '.csv'), delimiter=',')
        
        if max(data.flatten()) - min(data.flatten()) < 1e-6:
            if np.abs(np.mean(data) - 1) < 0.1:
                category = '1'
            else:
                category = '2'
            
        else:
            subprocess.run([filename], shell=True)
            category = input(USER_MENU)
        
        with open('classification.csv', 'a') as f:
            f.write(','.join((phi, gamma, category)) + '\n')
            
def plot_grid():
    df = pd.read_csv('classification.csv')
    df['category'] = df['category'].astype(int)
    plt.scatter(df['phi'], df['gamma'], c=df['category'], cmap='brg')
    plt.colorbar()
    plt.xlabel('phi')
    plt.ylabel('gamma')
    plt.show()
    
if __name__ == '__main__':
    # fill_grid()
    plot_grid()