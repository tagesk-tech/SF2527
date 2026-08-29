# SF2527

MATLAB solutions and LaTeX reports for SF2527 Numerical Methods for
Differential Equations I.

## Computer Exercise 1

- `Computer_exercise1/matlab/`: runnable MATLAB scripts and helper functions
- `Computer_exercise1/report/`: LaTeX source, figures, and compiled report
- `Computer_exercise1/CE1ht26.pdf`: assignment

Run a MATLAB part from `Computer_exercise1/matlab/`, for example:

```matlab
ce1_part1_rk_landau_lifshitz
```

Build the report with:

```bash
cd Computer_exercise1/report
cp student_info.example.tex student_info.tex
# Add the group details to student_info.tex. This file is ignored by Git.
make pdf
```

## Working together

```bash
git clone https://github.com/tagesk-tech/SF2527.git
cd SF2527
git switch -c name/topic
git pull --rebase origin main
git push -u origin name/topic
```

Open a pull request on GitHub when a change is ready to merge. Because the
repository is public, keep personal numbers and email addresses only in the
ignored `student_info.tex` file.
