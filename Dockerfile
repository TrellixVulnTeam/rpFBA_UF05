FROM brsynth/rpbase:dev

RUN pip install --no-cache-dir cobra
COPY rpTool.py /home/
COPY rpToolServe.py /home/
COPY galaxy/tool_rpFBA.py /home/
