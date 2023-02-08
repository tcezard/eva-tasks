FROM python:3.10

ENV vcf_validator_version=0.9.4
ENV python_version=3.10.9
ENV submission_automation_version=v1.9.0
ENV samtools_version=1.16.1
ENV bcfools_version=1.16
ENV htslib_version=1.16
ENV NXF_VER=22.10.6

WORKDIR /opt

RUN mkdir vcf_validator \
    && cd vcf_validator  \
    && curl -LJo vcf_validator  https://github.com/EBIvariation/vcf-validator/releases/download/v${vcf_validator_version}/vcf_validator_linux \
    && curl -LJo vcf_assembly_checker  https://github.com/EBIvariation/vcf-validator/releases/download/v${vcf_validator_version}/vcf_assembly_checker_linux \
    && chmod 755 vcf_assembly_checker vcf_validator \
    && cd ..

# Install python
#RUN curl -LJO https://www.python.org/ftp/python/${python_version}/Python-${python_version}.tgz
#RUN tar -zxf Python-${python_version}.tgz
#RUN rm Python-${python_version}.tgz
#RUN cd Python-${python_version}  \
#    && ./configure --prefix=$(pwd) > configure.log
#RUN cd Python-3.10.9 && make > make.log 2>&1
#RUN cd Python-3.10.9 && make altinstall > install.log 2>&1
#RUN cd ..
#RUN Python-${python_version}/python -m venv ./venv

RUN curl -LJO https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 \
    && tar xvf samtools-${samtools_version}.tar.bz2 \
    && cd samtools-${samtools_version} \
    && make install

RUN curl -LJO https://github.com/samtools/bcftools/releases/download/${bcfools_version}/bcftools-${bcfools_version}.tar.bz2 \
    && tar xvf bcftools-${bcfools_version}.tar.bz2 \
    && cd bcftools-${bcfools_version} \
    && make install

RUN curl -LJO https://github.com/samtools/htslib/releases/download/${htslib_version}/htslib-${htslib_version}.tar.bz2 \
    && tar xvf htslib-${htslib_version}.tar.bz2 \
    && cd htslib-${htslib_version} \
    && make install


RUN apt update && apt install -y default-jdk

RUN curl -s https://get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/


COPY config_file.yml /root/.submission_config.yml

RUN pip install -q "git+https://github.com/EBIvariation/eva-submission.git@${submission_automation_version}"


ENTRYPOINT ["validate_submission.py --eload 1"]
