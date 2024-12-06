#include <iostream>
#include <sundials/sundials_types.h>
#include <ida/ida.h>                // For the IDA solver
#include <nvector/nvector_serial.h>  // Serial N_Vector
#include <sunlinsol/sunlinsol_dense.h>  // Dense linear solver
// #include <ida/ida_dense.h>           // Dense solver interface
#include <sundials/sundials_math.h>  // SUNDIALS math utilities

// Function to define the system of equations (DAE)
int residuals(sunrealtype t, N_Vector y, N_Vector yp, N_Vector res, void *user_data) {
    sunrealtype y1 = NV_Ith_S(y, 0);  // y1
    sunrealtype y2 = NV_Ith_S(y, 1);  // y2
    sunrealtype yp1 = NV_Ith_S(yp, 0); // y1'

    // Differential equation: y1' = y2
    NV_Ith_S(res, 0) = yp1 - y2;

    // Algebraic equation: y1^2 + y2^2 = 1
    NV_Ith_S(res, 1) = y1 * y1 + y2 * y2 - 1;

    return 0;
}

int main() {
    // Problem size
    const int NEQ = 2;

    // Create the serial vectors for solution and derivatives
    N_Vector y = N_VNew_Serial(NEQ);   // Solution vector
    N_Vector yp = N_VNew_Serial(NEQ);  // Derivative vector
    N_Vector y0 = N_VNew_Serial(NEQ);  // Initial guess for y

    // Initialize the system (starting at point on the unit circle)
    NV_Ith_S(y, 0) = 1.0;  // y1 = 1
    NV_Ith_S(y, 1) = 0.0;  // y2 = 0

    NV_Ith_S(yp, 0) = 0.0;  // y1' = 0 (initial derivative guess)
    NV_Ith_S(yp, 1) = 1.0;  // y2' = 1 (since y1' = y2)

    // Create the IDA solver object
    void *ida_mem = IDACreate();
    if (ida_mem == NULL) {
        std::cerr << "Error: Unable to create IDA solver memory." << std::endl;
        return -1;
    }

    // Initialize the IDA solver
    sunrealtype t0 = 0.0;
    if (IDAInit(ida_mem, residuals, t0, y, yp) != IDA_SUCCESS) {
        std::cerr << "Error: Unable to initialize IDA solver." << std::endl;
        return -1;
    }

    // Set solver options
    IDASStolerances(ida_mem, 1.0e-6, 1.0e-8);

    // Create dense linear solver
    SUNLinearSolver LS = SUNNonlinearSolver(y, SUNLinSol_Dense(y));
    IDASetLinearSolver(ida_mem, LS, SUNLinSol_Dense(y));

    // Time-stepping
    sunrealtype t1 = 1.0;
    sunrealtype tret;
    while (t0 < t1) {
        int flag = IDASolve(ida_mem, t1, &tret, y, yp, IDA_NORMAL);
        if (flag < 0) {
            std::cerr << "Error during integration" << std::endl;
            break;
        }
        // Output the solution at each time step
        std::cout << "t = " << tret << ", y1 = " << NV_Ith_S(y, 0)
                  << ", y2 = " << NV_Ith_S(y, 1) << std::endl;
        t0 = tret;
    }

    // Free memory
    N_VDestroy(y);
    N_VDestroy(yp);
    IDAFree(&ida_mem);
    SUNLinSolFree(LS);

    return 0;
}
