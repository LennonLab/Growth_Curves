Bacterial Growth Curves
=============

## Lennon Lab Growth Curve Analysis Scripts

The code contained in the repo provides an anlaysis pipeline for bacterial growth curves.
The code fits a Modified Gompertz Equation to optical density readings.
The code uses grid searching maximum likelihood estimation to fit the model.

## How to use this code

1. Copy the contents of `bin` into your project
2. Use `test.example.R` as a template for your analysis
3. Proceed with analysis

## Repo Contents

**bin**

* *modified_Gomp.R*: function for analysis pipeline
* *read.synergy.R*: contains code for parsing raw data from synergy MX
* *grid.mle2.R*: contains code to grid search maximum likelihood
* *curve_fit_fxs.R*: contains the growth model equations		

**data**

* Various raw data examples

**test**

* *test.example.R*: A script with example code


## Contributors

[Mario Muscarella](http://mmuscarella.github.io/): Ph.D. candidate in the [Lennon Lab](http://www.indiana.edu/~microbes/people.php).

[Dr. Jay Lennon](http://www.indiana.edu/~microbes/people.php): Principle Investigator, Associate Professor, Department of Biology, Indiana University, Bloomington. Head of the [Lennon Lab](http://www.indiana.edu/~microbes/people.php).

Venus Kuo: Ph.D. student in the [Lennon Lab](http://www.indiana.edu/~microbes/people.php).

## License

Please see LICENSE.md
