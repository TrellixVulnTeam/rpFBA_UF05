FROM brsynth/rpcache:dev

RUN pip install --no-cache-dir cobra==0.16

RUN git clone https://github.com/Galaxy-SynBioCAD/inchikeyMIRIAM.git -b standalone
RUN mv inchikeyMIRIAM/inchikeyMIRIAM.py /home/
RUN rm -r inchikeyMIRIAM

COPY rpTool.py /home/
COPY rpToolServe.py /home/
COPY galaxy/code/tool_rpFBA.py /home/
