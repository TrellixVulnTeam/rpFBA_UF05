FROM ibisba/rpsbml

RUN pip install --no-cache-dir cobra

COPY rpFBA.py /home/
