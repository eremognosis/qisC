# qisC Error Codes

This file lists the error codes currently emitted by qisC. The source of truth is
the `setlasterror()` calls in `src/`.
The current implementation stores the last error in a thread-local 512-byte
buffer and prints each set error to `stderr` with a debug prefix. I thought of adding a "turn off button" to yelling errors but for now too hungry.
I wish there existed some auto doc maker in my flavour

## Families

| Prefix | Area |
| --- | --- |
| `ERRQBT` | Circuit and state-vector lifecycle |
| `ERRGT` | Gate creation, validation, and application |
| `ERRTNSR` | Tensor and state-vector math helpers |
| `ERRGRV` | Grover search helpers |

## Circuit Errors

| Code | Current message | Trigger | Recommended action |
| --- | --- | --- | --- |
| `ERRQBT001` | `n must be in 1 to 30, got %d` | `createcirc()` received a qubit count outside the supported range. | Pass a qubit count from `1` to `30`. Keep practical memory limits in mind. |
| `ERRQBT002` | `Memory Allocation Failed` | `createcirc()` could not allocate the `circuitt` struct. | Check available memory. |
| `ERRQBT003` | `Memory Allocation Failed for the state amplitudes` | `createcirc()` could not allocate the dense state vector. | Reduce qubit count or free memory. |
| `ERRQBT004` | `RESET CIRCUIT FAILED. NULL CIRCUIT` | `resetcirc()` received `NULL` or a circuit with no state buffer. | Pass a valid circuit created by `createcirc()`. |
| `ERRQBT005` | I forgot| `statevec()` received a `NULL` circuit, `NULL` state buffer, or `NULL` output buffer. | Pass a valid circuit and output buffer. |
| `ERRQBT006` | `STATEVECTOR EXPOPT BUFFER TOO SMALL. NEEDED %zu, got %zu` | `statevec()` output buffer length is smaller than `circuit->dim`. | Allocate at least `circuit->dim` `qScT` entries. |

## Gate Errors

| Code | Current message | Trigger | Recommended action |
| --- | --- | --- | --- |
| `ERRGT001` | `GATE WITH 0 DIMENSION` | `createGate(0)` was called. | Create gates with at least one qubit. |
| `ERRGT002` | `MEMORY ALLOCATION FAILED` | `createGate()` could not allocate the `Gate` struct. | Check available memory. |
| `ERRGT003` | `DATA MEMORY ALLOCATION FAILED` | Gate matrix allocation failed. | Reduce gate size or free memory. |
| `ERRGT004` | `NULL GATE` | A gate operation received `NULL`. | Pass a valid `Gate *`. |
| `ERRGT005` | `NULL DATA` or `NULL CIRCUIT` | `setGate()` received `NULL` matrix data, or `applygate()` received a `NULL` circuit. | Check the function context and pass the missing pointer. |
| `ERRGT006` | `NULL TARGET QUBITS` | `applygate()` received `NULL` target-qubit array. | Pass an array with one entry per gate qubit. |
| `ERRGT007` | `NULL GATE DATA` | A gate has no matrix data. | Use a gate from `createGate()` or a `CREATE_*()` factory and do not clear `gate->data`. |
| `ERRGT008` | `NULL CIRCUIT STATE` | `applygate()` received a circuit without a state buffer. | Use a valid circuit from `createcirc()`. |
| `ERRGT009` | `INVALID GATE QUBIT COUNT` | Gate qubit count is invalid or too large for `size_t` shifts. | Keep `gate->n` valid and consistent with `gate->dim`. |
| `ERRGT010` | `INVALID CIRCUIT QUBIT COUNT` | Circuit qubit count is invalid or too large for `size_t` shifts. | Use a valid circuit created by `createcirc()`. |
| `ERRGT011` | `GATE LARGER THAN CIRCUIT` | `gate->n > circuit->n`. | Apply only gates that fit inside the target circuit. |
| `ERRGT012` | `INVALID GATE DIMENSION` | `gate->dim` does not equal `2^gate->n`. | Do not manually corrupt `Gate` fields. |
| `ERRGT013` | `INVALID CIRCUIT DIMENSION` | `circuit->dim` does not equal `2^circuit->n`. | Do not manually corrupt `circuitt` fields. |
| `ERRGT014` | `TARGET QUBIT OUT OF RANGE` | A target index is negative or greater than/equal to `circuit->n`. | Keep target qubits in `[0, circuit->n - 1]`. |
| `ERRGT015` | `DUPLICATE TARGET QUBIT` | The same circuit qubit appears more than once in `target_qubits`. | Use unique target qubit indexes. |
| `ERRGT016` | `STATE MEMORY ALLOCATION FAILED` | Temporary allocation failed during gate application or unitary conversion. | Reduce circuit/gate size or free memory. |
| `ERRGT017` | `GATE MATRIX TOO LARGE` | Gate matrix size overflowed `size_t` or allocation-size calculations. | Use a smaller gate. |
| `ERRGT018` | `INVALID UNITARY CHECK EPSILON` | `checkunitary()` received a negative or non-finite epsilon. | Pass a finite epsilon greater than or equal to zero. |
| `ERRGT019` | `NON-FINITE GATE MATRIX VALUE` | Gate data contains `NaN` or infinity, or unitary checking produced a non-finite value. | Validate matrix data before calling `setGate()` or `makeunitary()`. |
| `ERRGT020` | `COULD NOT MAKE GATE UNITARY` | `makeunitary()` could not orthonormalize the matrix. | Provide a better-conditioned matrix or create a known unitary gate. |
| `ERRGT021` | `GATE IS NOT UNITARY` | `setGate()` rejected non-unitary matrix data. | Use a unitary matrix or call `makeunitary()` on a separate gate before setting. |

## Tensor Errors

| Code | Current message | Trigger | Recommended action |
| --- | --- | --- | --- |
| `ERRTNSR001` | `STATEVECTOR IS NULL` | Tensor helper received a `NULL` circuit pointer. | Pass a valid circuit. |
| `ERRTNSR002` | `STATEVECTOR IS EMPTY` | Tensor helper received a circuit whose `state_amps` is `NULL`. | Use a valid circuit from `createcirc()`. |
| `ERRTNSR003` | `STATE VECTOR OF DIFFERENT DIMENSIONS` | `tensor_inner_product()` received circuits with different dimensions. | Compare states with equal `dim`. |
| `ERRTNSR004` | `INVALID INPUT` | Tensor helper received invalid matrix pointers, output pointers, or dimensions. | Check all pointers and dimensions before calling. |
| `ERRTNSR006` | `ZERO VECTOR` | `tensor_normalize()` cannot normalize a zero vector. | Ensure at least one amplitude is non-zero before normalizing. |

(005 is currently unused mostly because highly caffeinaated coders at 3 am is bad at counting.)

## Grover Errors

| Code | Current message | Trigger | Recommended action |
| --- | --- | --- | --- |
| `ERRGRV001` | `NULL CIRCUIT` | Grover helper received `NULL` circuit. | Pass a valid circuit. |
| `ERRGRV002` | `NULL CIRCUIT STATE` | Circuit has no state buffer. | Use a valid circuit from `createcirc()`. |
| `ERRGRV003` | `INVALID CIRCUIT QUBIT COUNT` | Circuit qubit count is invalid or too large for `size_t` shifts. | Use a valid circuit and avoid manual field edits. |
| `ERRGRV004` | `INVALID CIRCUIT DIMENSION` | `circuit->dim` does not equal `2^circuit->n`. | Do not manually corrupt `circuitt` fields. |
| `ERRGRV005` | `NULL ORACLE` | Oracle callback is `NULL`. | Pass a valid `grover_mark_fn`. |
| `ERRGRV006` | `ORACLE MARKED NO STATES` | Oracle did not mark any basis state. | Fix the oracle or target state. |
| `ERRGRV007` | `NON-FINITE STATE AMPLITUDE` | State validation found `NaN` or infinity. | Check previous operations and input state. |
| `ERRGRV008` | `TARGET MEMORY ALLOCATION FAILED` | Could not allocate internal oracle marked-state cache. | Reduce qubit count or free memory. |
| `ERRGRV009` | `INVALID QUBIT COUNT` | `groverrecit()` received an invalid qubit count. | Pass a positive qubit count that fits in `size_t` shifts. |
| `ERRGRV010` | `INVALID MARKED STATE COUNT` | Marked count is zero or larger than the search space. | Pass `1 <= marked_count <= 2^n`. |
| `ERRGRV011` | `RECOMMENDED ITERATION COUNT TOO LARGE` | Estimated Grover iteration count exceeds `INT_MAX`. | Use a smaller problem or choose iterations manually. |
| `ERRGRV012` | `MARKED STATE OUT OF RANGE` | `grover_search()` marked state is greater than/equal to `circuit->dim`. | Use a basis state in `[0, circuit->dim - 1]`. |
| `ERRGRV013` | `NULL OUTPUT STATE` | `grovermlstate()` received `NULL` for `out_state`. | Pass a valid `size_t *`. |

## Notes For Maintainers(Me)
-- no notes im good
-- there will be more errors as we (I) add features anyway so lemme pet the cat and study chem
