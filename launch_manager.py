import subprocess
import time


def arange(start, stop, step):
    while start < stop:
        yield start
        start += step
        
MAX_SUBMITS = 400

U_s         = 0.4
phi_tab     = list(arange(0.5, 1.0, 0.01))
gamma_tab   = list(arange(0.2, 0.55, 0.01))

params = {(phi, gamma) for phi in phi_tab for gamma in gamma_tab}

while len(params) > 0:
    
    to_launch = min(len(params), MAX_SUBMITS-(subprocess.run(['squeue', '--me'], stdout=subprocess.PIPE).stdout.decode("utf-8").count("\n")-1))
    
    for i in range(to_launch):
        phi, gamma = params.pop()
        print("Launching job with U_s = ", U_s, ", phi = ", phi, " and gamma = ", gamma)
        subprocess.run(['sbatch', 'single_test.sh',
                        '-U', str(U_s),
                        '-p', str(phi),
                        '-g', str(gamma)
        ])
        time.sleep(0.01)

    time.sleep(10)
    
print("All jobs launched")