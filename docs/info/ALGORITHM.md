# Digital Multi-Load Beam Structural Analysis and Simulation System

## Program algorithm

**Step 1.** Start the program and display the beam input form.

**Step 2.** Enter the beam length. Choose **Constant EI** and enter EI, or choose **E and I** and enter the material stiffness E and section property I. Both choices use one constant stiffness along the beam.

**Step 3.** Choose support types and their positions. Add the point loads, uniformly distributed loads, triangular or trapezoidal loads, and applied moments. Select a row to edit or remove it when needed.

**Step 4.** Click **Calculate**. Check that all active inputs are finite, the length and stiffness are positive, supports restrain the beam, and all supports and loads lie within the beam. If an input is invalid, display an error and return to the input form.

**Step 5.** Calculate the support reactions, shear force, bending moment, slope and deflection. Check that the calculation is reliable. If it cannot be solved reliably, display an error and let the user correct the model.

**Step 6.** Display the reactions, shaded graphs and magnified deformed beam shape. Show the largest results and their positions.

**Step 7.** Enter a position x and click **Read** to see the shear, moment, slope and deflection there. Click **Show Solution** to view the formulas and calculation details.

**Step 8.** If the model is changed, clear the previous results and return to the input steps. Otherwise, close the program to end.
