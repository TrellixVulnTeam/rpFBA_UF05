FROM brsynth/rpbase

RUN pip install --no-cache-dir cobra && \
    conda install -c conda-forge flask-restful

COPY rpFBA.py /home/
COPY rpFBAServe.py /home/

ENTRYPOINT ["python"]
CMD ["/home/rpFBAServe.py"]

# Open server port
EXPOSE 8994
