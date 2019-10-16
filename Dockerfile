FROM brsynth/rpbase

WORKDIR /home/

RUN pip install --no-cache-dir cobra

COPY rpFBA.py /home/
