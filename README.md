# japan_ferns_review

Code repository for Ebihara and Nitta. "An update and reassessment of fern and lycophyte diversity data in the Japanese Archipelago", Journal of Plant Research (2019). https://doi.org/10.1007/s10265-019-01137-3

All code is in [R](https://cran.r-project.org/). The [drake package](https://ropensci.github.io/drake/) is used to manage the workflow.

Raw data (except for a single file, `FernGreenListV1.01.xls`, which is included in this repo) are kept on [Dryad](https://datadryad.org/) (https://doi.org/10.5061/dryad.4362p32).

To run all analyses and generate the manuscript, first [clone this repository](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository). Next, [download the raw data from the Dryad repository](https://doi.org/10.5061/dryad.4362p32), unzip it, and move it to the `data/` folder. Finally, run `make.R`.

## Reproducible analysis with Docker

`make.R` requires various packages to be installed, and may not work properly if package versions have changed. Therefore, a [Docker image is provided](https://hub.docker.com/r/joelnitta/japan_ferns_review) to run the code reproducibly.

To use it, first [install docker](https://docs.docker.com/install/).

Clone this repo and download the data as described above.

Navigate to the cloned repository (where `/path/to/repo` is the path on your machine), and launch the container:

```
cd /path/to/repo
docker-compose up -d
```

If `docker-compose` is not available, you can also launch the container using `docker run`. Be sure to replace the path on the left side of `:` (in this example, `/path/to/japan_ferns_review/`) with the path to this repo on your machine.

```
docker run -d --name japan_ferns_review_analysis_1 -e DISABLE_AUTH=true -v /path/to/japan_ferns_review/:/home/rstudio/japan_ferns_review joelnitta/japan_ferns_review:3.6.0 bash
```

Enter the container:

```
docker exec -it japan_ferns_review_analysis_1 bash
```

Inside the container, run `make.R`:

```
Rscript make.R
```

You will see the targets being built by `drake`, and the final manuscript and figures should be compiled at the end as in the `manuscript` folder as `japan_ferns_review_ms.docx`, `fig_1.pdf`, etc.

When it's finished, exit the container and take it down:

```
exit
docker-compose down
```
