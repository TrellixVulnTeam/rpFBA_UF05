FROM brsynth/rpbase

WORKDIR /home/

RUN pip install --no-cache-dir cobra && \
    conda install -c conda-forge flask-restful

COPY rpFBA.py /home/
COPY rpFBAServe.py /home/

ENTRYPOINT ["python"]
CMD ["/home/rpFBAServe.py"]

EXPOSE 5005
