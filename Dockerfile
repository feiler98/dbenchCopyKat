FROM feiler98/copykat_pyomics

RUN mkdir -p /scratch/tmp/feiler/dbenchCopyKat
WORKDIR /scratch/tmp/feiler/dbenchCopyKat
COPY . .

RUN yum install -y pip
RUN yum install -y python3-devel
RUN yum install -y libpng-devel
RUN pip install --no-cache-dir -r requirements.txt
RUN R -e "install.packages('Seurat',dependencies=TRUE, repos='http://cran.rstudio.com/')"

CMD ["python3", "/scratch/tmp/feiler/dbenchCopyKat/run_copykat.py"]