gunicorn app:server --chdir /primeDesign/web_app --bind 0.0.0.0:9994 --timeout 1800 --access-logfile - --workers 4
# python /primeDesign/web_app/app.py