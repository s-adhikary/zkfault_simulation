To build the benchmarking binaries type `./build.sh`.

To clean the built files type `./clean`

Go to `bin` directory to find the executatble binaries. To change the simulation parameters go to the file `Test_codes/lib/less_test_fault.c` and set the parameters as macro. Here are the following parameters:

- `PROB`: The probability of successful fault injection. Set it as a `float` or `double`.
- `CASE`: There are three cases:
    - `CASE 1`: Run simulation to count the recovered secret columns with single fault.
    - `CASE 2`: Run simulation to count the faulted signatures to recover complete secret.
    - `CASE 3`: Same as `CASE 2` but it also takes `PROB` into account.
- `NUM_TEST_ITERATIONS`: Total number of time we run the simulation.

Change the location of fault in `Reference_Implementation/lib/seedtree.c` file.
