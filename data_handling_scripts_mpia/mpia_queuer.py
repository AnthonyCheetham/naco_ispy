"""

This script is intended to run in the background. It reads in a list of jobs 
to be executed from a text file and launches them in the background when the 
cluster is not already being overused.

It should work in tandem with the naco_ispy queue_analysis.py version made 
for MPIA. 

You can run a screen instance with 
    screen -R ISPYQueue
then launch it with
    python mpia_queuer.py
and disconnect the screen window with ctrl+a, d
then it should run in the background. You can check on it with
    screen -R ISPYQueue

"""
import psutil
import numpy as np
import time
import subprocess
import datetime
import os


check_interval = 60 # sec
queue_file = '/data/beegfs/astro-storage/groups/henning/rlau/naco-ispy-shared/ispy_job_queue.que'
idle_cpu_percent = 50.0 # if it is less than this, the script will launch a job

# Read/write functions for the queue file
def read_file(queue_file):
    jobs = np.loadtxt(queue_file,delimiter=",",dtype=str)
    return jobs

def write_file(queue_file,jobs):
    np.savetxt(queue_file,jobs,delimiter=",",fmt="%s")

# 
def main_loop():
	print('ISPY job queuer started!')
	print('Press ctrl+c to exit')
	while True:
		try:
			# Check queue file for jobs in the queue
			jobs = read_file(queue_file)
            if jobs.ndim == 1:
                jobs = np.array([jobs])

			# Check cpu usage over 10 secs
			cpu_usage = psutil.cpu_percent(interval=10)

			# 
			if (cpu_usage < idle_cpu_percent) and (jobs.size > 0):
				next_job = jobs[0]

				# change into the directory and launch the job in the background
				directory,next_job_cmd = next_job
				os.chdir(directory)

				cmd = 'echo "'+next_job_cmd+'" | at now' # i.e. launch in background

				# Then launch the next job
				print('{0}: Launching: {1}'.format(datetime.datetime.now(),next_job))
				subprocess.call(cmd,shell=True)

				# Then save the updated job list
				write_file(queue_file,jobs[1:])

			# Some outputs for debugging
			# elif (len(jobs) == 0):
				# print('Nothing to do')
			# else:
				# print('CPU usage too high: {0}'.format(cpu_usage))

			# Sleep 
			time.sleep(check_interval)

		except KeyboardInterrupt:
			print('Exiting ISPY job queuer')
			break
# This is so it runs if you type python mpia_queuer.py, but can also be imported as a module 
if __name__ == "__main__":
	# Run the loop
	main_loop()