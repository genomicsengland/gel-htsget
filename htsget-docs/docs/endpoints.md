# Endpoints

## Overview

* **Login** - `POST /user/login`
    - Given a LDAP username and a password (`{"username": "USER", "password": "PASS"}`), returns a session ID.
    - This token is passed in the header of the GET requests to access the genomic data.
    - Authentication and authorisation are managed by [OpenCGA](http://docs.opencb.org/display/opencga).

---

* **Reads** - `GET /reads/{study_id}/{sample_id}`
    - This endpoint allows access to read data (BAM)

---

* **Variants** - `GET /variants/{study_id}/{vcf_type}/{sample_id}`
    - This endpoint allows access to variant data (VCF)
    - The URL parameter `vcf_type` refers to the type of data stored in the VCF (STR, CNV, SV, small_variant, family, somatic) 


## Query parameters

To know more about the parameters that can be used to filter the data, please refer to [htsget specification](http://samtools.github.io/hts-specs/htsget.html).