# japan_ferns_review

Code repository and raw data for Ebihara and Nitta. "An update and reassessment of fern and lycophyte diversity data in the Japanese Archipelago", Journal of Plant Research.

All code is in [R](https://cran.r-project.org/). The [drake package](https://ropensci.github.io/drake/) is used to manage the workflow. To run all analyses and generate the manuscript, simply [clone this repository](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository) and run `make.R`.

## Reproducible analysis with Docker

`make.R` requires various packages to be installed, and may not work properly if package versions have changed. Therefore, a [Docker image is provided](https://hub.docker.com/r/joelnitta/japan_ferns_review) to run the code reproducibly.

To use it, first [install docker](https://docs.docker.com/install/) and clone this repository.

Navigate to the cloned repository (where `/path/to/repo` is the path on your machine), and launch the container:

```
cd /path/to/repo
docker-compose up -d
```

Enter the container:

```
docker exec -it japan_ferns_review_analysis_1 bash
```

Inside the container, run `make.R`:

```
Rscript make.R
```

You will see the targets being built by `drake`, and the final manuscript and figure should be compiled at the end as in the `manuscript` folder `japan_ferns_diversity_ms.docx` and `fig_1.pdf`, respectively.

When it's finished, exit the container and take it down:

```
exit
docker-compose down
```