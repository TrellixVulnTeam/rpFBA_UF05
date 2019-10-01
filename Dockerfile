FROM brsynth/rpbase

RUN pip install --no-cache-dir cobra

COPY rpFBA.py /home/
