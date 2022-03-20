#! /usr/bin/python


import subprocess, threading


class Command(object):
    def __init__(self, cmd):
        self.cmd = cmd
        self.process = None

    def run(self, timeout):
        def target():
            print 'Thread started'
            self.process = subprocess.Popen(self.cmd, shell=True)
            self.process.communicate()
            print 'Thread finished'

        thread = threading.Thread(target=target)
        thread.start()

        thread.join(timeout)
        if thread.is_alive():
            print 'Terminating process'
            self.process.terminate()
            thread.join()

	returncode = self.process.returncode
        print returncode

	return returncode


if __name__ == "__main__":

	command = Command("/home/asassn/software/diapl/fwhm_ping/fwhm fwhm.par timeout_example_input.fits") 
	command.run(timeout=50)
