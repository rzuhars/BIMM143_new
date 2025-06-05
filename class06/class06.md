# Class 06: R Functions
Renee Zuhars (PID: A17329856)

- [Section 1. Writing Functions](#section-1-writing-functions)
  - [Every R function has 3 things:](#every-r-function-has-3-things)
  - [Using the sample() function](#using-the-sample-function)
- [Section 2. Generate DNA sequence](#section-2-generate-dna-sequence)
  - [Using the `together=` logical](#using-the-together-logical)
- [Section 3. Generate Protein
  Function](#section-3-generate-protein-function)
  - [I found this solution using
    chatgpt:](#i-found-this-solution-using-chatgpt)
  - [In class, we used this solution (using
    sapply()):](#in-class-we-used-this-solution-using-sapply)

# Section 1. Writing Functions

Let’s start writing our first silly function to add some numbers.

### Every R function has 3 things:

- name (we get to pick this)
- input arguments (there can be loads of these separated by a comma)
- the body (the R code that does the work)

> Note: before modifications later in this exercise, this function read:
> add \<- function(x,y){x + y}

``` r
add <- function(x, y=100, z=0){
  x + y + z
}
```

I can just use this function like any other function as long as R knows
about it- which means I have to make sure to run the previous code chunk
first:

``` r
add(1, 100)
```

    [1] 101

``` r
add(x=c(1,2,3,4), y=100)
```

    [1] 101 102 103 104

What if we only put one variable inside the function?

> In order to do this, we need to add a default to the original
> function. The original function has been modified to set a default of
> y=100.

``` r
add(1)
```

    [1] 101

Functions can have “required” input arguments and “optional” input
arguments. The optional arguments are defined with an equals default
value (`y=100`) in the function definition.

> Here we have added another variable. The original function was
> modified again to set a default of z=0. Notice how despite this
> modification, the other code chunks still work!

``` r
add(x=1, y=100, z=10)
```

    [1] 111

### Using the sample() function

> Q. Write a function to return a DNA sequence of a user specified
> length. Call it `generate_dna()`

The `sample()` function can help here:

``` r
#generate_dna <- function(size=5){}

students <- c("jeff","jeremy","peter")

sample(students, size=5, replace=TRUE)
```

    [1] "jeff"   "jeremy" "jeremy" "peter"  "jeremy"

Above, the “replace” argument was used to avoid the error of asking for
more than the population size.

# Section 2. Generate DNA sequence

> Now work with `bases` rather than `students`

``` r
bases <- c("A", "C", "G", "T")

sample(bases, size=10, replace=TRUE)
```

     [1] "A" "T" "C" "A" "T" "C" "T" "T" "A" "T"

Now that I have a working ‘snippet’ of code, I can use this as the body
of my first function version here.

> Below, I have changed the ‘size’ parameter so that it is not limited
> to only 10 characters. Now I can generate a 100 character long
> sequence (the size=5 portion of the code is the default, and make the
> default negligible by adding size=size):

``` r
generate_dna <- function(size=5) {
  bases <- c("A", "C", "G", "T")
  sample(bases, size=size, replace=TRUE)
}

generate_dna(100)
```

      [1] "T" "T" "A" "A" "C" "T" "C" "A" "T" "T" "G" "T" "C" "T" "C" "T" "C" "T"
     [19] "C" "G" "T" "T" "A" "T" "A" "A" "G" "G" "C" "T" "G" "A" "A" "T" "G" "C"
     [37] "T" "G" "G" "G" "C" "A" "T" "A" "G" "T" "C" "A" "G" "A" "G" "G" "A" "C"
     [55] "A" "T" "C" "T" "A" "T" "T" "C" "C" "C" "G" "A" "C" "C" "C" "T" "A" "C"
     [73] "C" "A" "C" "T" "C" "C" "T" "A" "C" "T" "T" "A" "A" "G" "T" "C" "T" "C"
     [91] "G" "A" "A" "G" "A" "G" "C" "T" "C" "G"

### Using the `together=` logical

I want the ability to return a sequence like “AGTACCTG” - I want all the
characters in one element vector where the bases are all together (not
spaced out like above), so it can be pasted into BLAST or some other
tool…

``` r
generate_dna <- function(size=5, together=TRUE) {
  bases <- c("A", "C", "G", "T")
  sequence <- sample(bases, size=size, replace=TRUE)
  
  if(together) {
    sequence <- paste(sequence, collapse="")
  }
  return(sequence)
}

generate_dna()
```

    [1] "AAGGG"

The sequence above is generated in the format that we wanted. To undo
this, we can use the “together” logical:

``` r
generate_dna(together=FALSE)
```

    [1] "T" "T" "G" "G" "A"

# Section 3. Generate Protein Function

> Q. Write a function, `generate_protein()`, to return protein sequences
> of user determined length.

We can get the set of 20 natural amino acids from the **bio3d** package.

``` r
aa <- bio3d::aa.table$aa1[1:20]
```

The \$aa1 only returns the one-letter code of the 20 natural AAs.

``` r
generate_protein <- function(size=20, together=TRUE) {
  ## Get the 20 amino acids as a vector
   amino_acids <- c(aa)
  sequence <- sample(amino_acids, size=size, replace=TRUE)
  
  ## Optionally return a single element string
  if(together) {
    sequence <- paste(sequence, collapse="")
  }
  return(sequence)
}

generate_protein()
```

    [1] "FYRHAPEGGVIPTLVIWTSC"

> Q. Generate random protein sequences of length 6 to 12 amino acids

### I found this solution using chatgpt:

``` r
num_vars <- sample(6:12, size=1)

generate_protein <- function(size=num_vars, together=TRUE) {
  amino_acids <- c(aa)
  sequence <- sample(amino_acids, size=size, replace=TRUE)
  
  if(together) {
    sequence <- paste(sequence, collapse="")
  }
  return(sequence)
}

generate_protein()
```

    [1] "FRLPMYA"

The num_vars function selects a random integer between 6 and 12, then
uses that integer as the size, returning one sequence at a time of
random length btw 6-12 AAs.

### In class, we used this solution (using sapply()):

We can fix the inability to generate multiple sequences by either
editing and adding to the function body code (e.g. a for loop) or by
using the R **apply** family of utility function.

``` r
## Using the sapply function, we set the first argument equal to a vector consisting of the lengths we want to generate. The second argument is the function in question. 

sapply(6:12, generate_protein)
```

    [1] "FMVDGS"       "LKLAPFR"      "FLFWTVDT"     "NPEEAYMTD"    "PNPTDTWYYK"  
    [6] "DNPQPSAKCHM"  "DGPQPSPPRCYL"

It would be cool and useful if I could get FASTA format output.

``` r
ans <- sapply(6:12, generate_protein)
ans
```

    [1] "MRGNTC"       "TQLNWMD"      "HTKDWFNQ"     "AGTQTLYQL"    "SGLCARSDFP"  
    [6] "MQRSFGRYFVL"  "HETPRVWANWWA"

``` r
cat(ans, sep="\n")
```

    MRGNTC
    TQLNWMD
    HTKDWFNQ
    AGTQTLYQL
    SGLCARSDFP
    MQRSFGRYFVL
    HETPRVWANWWA

I want this to look like FASTA format with an ID line. The functions
`paste()` and `cat()` can help us here…

``` r
id.line <- paste(">ID.", 6:12, sep="")
seq.line <- paste(id.line, ans, sep="\n")
cat(seq.line, sep="\n")
```

    >ID.6
    MRGNTC
    >ID.7
    TQLNWMD
    >ID.8
    HTKDWFNQ
    >ID.9
    AGTQTLYQL
    >ID.10
    SGLCARSDFP
    >ID.11
    MQRSFGRYFVL
    >ID.12
    HETPRVWANWWA

> Q. Determine if these sequences can be found in nature or are they
> unique? Why or why not?

Simply copy/paste the above results into protein BLAST. (the results
used were from an iteration that resulted in the following data):

“\>ID.6 RFIDFA” “\>ID.7 TKRGLGF” “\>ID.8 RSQLENYY” “\>ID.9 THMMFWWQQ”
“\>ID.10 KENIMVNPQE” “\>ID.11 CKKGKKMDNDL” “\>ID.12 NLQALYNHWHAT”

I BLASTp searched my FASTA format sequences against refseq_protein, and
found that lengths 6, 7, and 8 are not unique and can be found in the
databases with 100% coverage and 100% identity. Random sequence lengths
9, 10, 11, and 12 are unique, and no matches had both 100% coverage and
100% identity within the database.
