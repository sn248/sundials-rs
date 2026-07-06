use thiserror::Error;

#[derive(Debug, Error)]
pub enum SundialsError {
    #[error("memory allocation failed in {0}")]
    Memory(&'static str),

    #[error("{solver} returned error code {code} in {function}")]
    ReturnCode {
        solver: &'static str,
        function: &'static str,
        code: i32,
    },

    #[error("right-hand side function returned non-zero flag {0}")]
    RhsError(i32),

    #[error("N_Vector creation failed")]
    NVectorCreate,

    #[error("SUNMatrix creation failed")]
    MatrixCreate,

    #[error("SUNLinearSolver creation failed")]
    LinearSolverCreate,

    #[error("step size became too small")]
    TooMuchWork,

    #[error("solution appears to diverge")]
    Diverge,
}

/// Map a SUNDIALS integer flag to `Ok(())` or an appropriate `Err`.
pub(crate) fn check_flag(
    flag: std::os::raw::c_int,
    solver: &'static str,
    function: &'static str,
) -> Result<(), SundialsError> {
    // Flag values are defined per-solver in their respective headers.
    // 0 is CV_SUCCESS / IDA_SUCCESS for all solvers.
    if flag >= 0 {
        Ok(())
    } else {
        Err(SundialsError::ReturnCode { solver, function, code: flag })
    }
}
