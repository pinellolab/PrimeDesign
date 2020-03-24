import os
import subprocess as sb
import sys

if len(sys.argv)==1:

	print('\n- PrimeDesign Docker Container -\n\t- use `primedesign_cli` for the command line tool\n\t- use `primedesign_webapp` for the web application\n')
	sys.exit(1)

if sys.argv[1]=='primedesign_cli':
	sb.call(['/opt/conda/bin/python', '/primeDesign/command_line/primeDesign.py']+ sys.argv[2:])
elif sys.argv[1]=='primedesign_webapp':
	sb.call(["/bin/bash", "-c", "/primeDesign/web_app/start_server_docker.sh"])
else:
	print('\n- PrimeDesign Docker Container -\n\t- use `primedesign_cli` for the command line tool\n\t- use `primedesign_webapp` for the web application\n')
