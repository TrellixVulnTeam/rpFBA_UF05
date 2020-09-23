FROM brsynth/rpcache:v2

#RUN pip install --no-cache-dir cobra==0.16
RUN rm -rf /usr/local/lib/python3.7/site-packages/ruamel*
RUN pip install --no-cache-dir cobra==0.16 timeout-decorator

RUN git clone https://github.com/Galaxy-SynBioCAD/inchikeyMIRIAM.git -b standalone
RUN mv inchikeyMIRIAM/inchikeyMIRIAM.py /home/
RUN rm -r inchikeyMIRIAM

COPY rpTool.py /home/
COPY rpToolServe.py /home/
COPY galaxy/code/tool_rpFBA.py /home/
