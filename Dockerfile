FROM brsynth/rpbase

RUN pip install --no-cache-dir cobra
COPY rpTool.py /home/
COPY rpToolServe.py /home/
COPY galaxy/code/tool_rpFBA.py /home/
