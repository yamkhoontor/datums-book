# Datums, Coordinate Systems and Map Projections

**A Practical Guide for Surveyors and Engineers, with SVY21 Worked Examples**
by Yam Khoon Tor

This repository holds the companion material of the book: the browser apps, the Python programs and the notebooks.

## The apps (no installation)

Open **https://yamkhoontor.github.io/datums-book/** in a web browser. There is one app for each chapter, with a button for each worked example in the chapter. The apps also run offline: press the green **Code** button, choose **Download ZIP**, unzip it and open `apps/index.html`.

## The notebooks in Google Colab (no installation)

Press a badge to open the chapter's notebook in Google Colab (a Google account is needed). Choose *Runtime → Run all*. Changes are not kept when the tab is closed, so the notebooks always reopen as published.

| Notebook | Open |
|---|---|
| Appendix E Test points | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YamKhoonTor/datums-book/blob/main/notebooks/Appendix_E_Test_points.ipynb) |
| Chapter 2 The Earth's figure | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YamKhoonTor/datums-book/blob/main/notebooks/Chapter_02_The_Earths_figure.ipynb) |
| Chapter 3 Coordinate systems | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YamKhoonTor/datums-book/blob/main/notebooks/Chapter_03_Coordinate_systems.ipynb) |
| Chapter 4 Heights and vertical datums | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YamKhoonTor/datums-book/blob/main/notebooks/Chapter_04_Heights_and_vertical_datums.ipynb) |
| Chapter 5 Classical datums and their realisation | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YamKhoonTor/datums-book/blob/main/notebooks/Chapter_05_Classical_datums_and_their_realisation.ipynb) |
| Chapter 7 ITRF and a moving Earth | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YamKhoonTor/datums-book/blob/main/notebooks/Chapter_07_ITRF_and_a_moving_Earth.ipynb) |
| Chapter 9 Transformation models | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YamKhoonTor/datums-book/blob/main/notebooks/Chapter_09_Transformation_models.ipynb) |
| Chapter 10 Deriving transformation parameters | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YamKhoonTor/datums-book/blob/main/notebooks/Chapter_10_Deriving_transformation_parameters.ipynb) |
| Chapter 11 Cross-nation datum transformation | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YamKhoonTor/datums-book/blob/main/notebooks/Chapter_11_Cross-nation_datum_transformation.ipynb) |
| Chapter 13 Transverse Mercator and the SVY21 grid | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YamKhoonTor/datums-book/blob/main/notebooks/Chapter_13_Transverse_Mercator_and_the_SVY21_grid.ipynb) |
| Chapter 14 Other projections a surveyor meets | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YamKhoonTor/datums-book/blob/main/notebooks/Chapter_14_Other_projections_a_surveyor_meets.ipynb) |
| Chapter 15 Grid to ground | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YamKhoonTor/datums-book/blob/main/notebooks/Chapter_15_Grid_to_ground.ipynb) |

Run the Appendix E notebook first on a new installation: it checks every function against published test points and reports PASS or FAIL.

## The programs

`programs/` has one Python program for each chapter with listings, and `appendix_e_tests.py`. Each is complete in itself and prints the output shown in the book. Python 3.9 or later; NumPy for Chapters 4, 9, 10 and 11 and for the tests.

```
python programs/ch13_svy21_tm.py
```

## Sources

Coastlines in the apps: OpenStreetMap data (© OpenStreetMap contributors, ODbL), via geo-maps. Plate boundaries in the book: Bird (2003), PB2002. Projection and transformation parameters: SLA, JUPEM, BIG, the Survey Department of Brunei and the EPSG registry, as cited in the book.

Copyright © 2026 Yam Khoon Tor. All rights reserved. See the copyright page of the book.
