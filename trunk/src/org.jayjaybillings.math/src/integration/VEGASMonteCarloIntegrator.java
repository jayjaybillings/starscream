/*******************************************************************************
 * Copyright (c) 2012-, Jay Jay Billings
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the author nor the
 *       names of the contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 * ARE DISCLAIMED. IN NO EVENT SHALL JAY JAY BILLINGS OR THE CONTRIBUTORS BE 
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 * POSSIBILITY OF SUCH DAMAGE.
 *******************************************************************************/

package integration;

import random.CMWCRNG;
import random.RNG;
import common.IFunction;

/** 
 * <!-- begin-UML-doc -->
 * <p>This class provides an implementation of the VEGAS algorithm invented by G. P. Lepage and is derived from the implementation in Numerical Recipes, 3rd Ed., Webnote No. 9, Rev. 1. Lepage's paper was used for reference: "A New Algorithm for Multidimensional Integration," SLAC-PUB-1839 (Revised) November 1976, revised April 1977.</p><p>This implementation has been modified from the original to make it into a class and break apart the original vegas() routine into smaller operations. The function pointer used in the original has been replaced by the IFunction interface.</p><p>There are four public, constant and static attributes of this class - NEW_GRID, NEW_ANSWERS, REUSE_ANSWERS and CONTINUE - the uses of which are explained in the documentation for setMode().</p><p>The VEGAS algorithm works by importance sampling and stratified sampling. It iteratively refines the grid for both methods up to a number of maximum iterations and making up to a certain number of calls to the function during each iteration. It converges on the answer much quicker than a simple Monte Carlo routine for most functions, but there are a few cases where it will not perform any better.</p>
 * <!-- end-UML-doc -->
 * @author jaybilly
 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
 */
public class VEGASMonteCarloIntegrator {
	/** 
	 * <!-- begin-UML-doc -->
	 * <p>This mode directs the integrator to construct a new grid and set all values to their defaults. This is the default mode and should be set with setMode() any time the integrator is going to be used to integrate a function from scratch.</p>
	 * <!-- end-UML-doc -->
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public static final int NEW_GRID = 0;
	/** 
	 * <!-- begin-UML-doc -->
	 * <p>This mode calculates new answers for the integration, but does so using a previously created grid. It resets all answers from previous runs, but uses the "optimized" grid. It can be very useful when an optimal grid for an integration has been calculate using the NEW_GRID mode and the integral needs to be evaluated to a higher accuracy. </p>
	 * <!-- end-UML-doc -->
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public static final int NEW_ANSWERS = 1;
	/** 
	 * <!-- begin-UML-doc -->
	 * <p>This mode reuses both the grid and the answers from calculations performed in other modes. However, it still resets the stratified sampling grid and rebins if needed. It is useful when integrating the same function multiple times, but the client wants the stratification reset.</p>
	 * <!-- end-UML-doc -->
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public static final int REUSE_ANSWERS = 2;
	/** 
	 * <!-- begin-UML-doc -->
	 * <p>This mode continues the integration from its previous state, changing none of the values from the last integration. It is useful for refining the grid for an additional n iterations without resetting the stratification.</p>
	 * <!-- end-UML-doc -->
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public static final int CONTINUE = 3;
	/** 
	 * <!-- begin-UML-doc -->
	 * <p>The operating mode in which the integration will be performed. By default it is equal to the NEW_GRID mode.</p>
	 * <!-- end-UML-doc -->
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	private int init_mode = 0;
	/** 
	 * <!-- begin-UML-doc -->
	 * <p>The maximum number of calls that will be made to the function during the integration. This value is one by default since there is no way to know how many calls a given client would prefer.</p>
	 * <!-- end-UML-doc -->
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	private int maxCalls = 1;
	/** 
	 * <!-- begin-UML-doc -->
	 * <p>The maximum number of iterations that the integrator should perform. By default this value is one since there is no </p>
	 * <!-- end-UML-doc -->
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	private int maxIterations = 1;
	/** 
	 * <!-- begin-UML-doc -->
	 * <p>An instance of an org.jayjaybillings.math.random.RNG class. By default, this class creates an instance of CMWCRNG.</p>
	 * <!-- end-UML-doc -->
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	private RNG randomNumber;
	/** 
	 * <!-- begin-UML-doc -->
	 * <p>The maximum number of divisions along a particular dimension. It is denoted K in Numerical Recipes and NDMX in the VEGAS code listing. It is set to 50 by default.</p>
	 * <!-- end-UML-doc -->
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	private final int maxAxialDivisions = 50;
	/** 
	 * <!-- begin-UML-doc -->
	 * <p>The maximum number of dimensions allowed in the integration. It is denoted d in Numerical Recipes and MXDIM in the VEGAS code listing. It is set to 10 by default.</p>
	 * <!-- end-UML-doc -->
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	private final int maxDimensions = 10;
	/** 
	 * <!-- begin-UML-doc -->
	 * <p>The convergence parameter. According to Lepage's paper, this value is typically set between 0.2 and 2.0 and controls the convergence of the grid.</p>
	 * <!-- end-UML-doc -->
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	private final double alpha = 1.5;
	/** 
	 * <!-- begin-UML-doc -->
	 * <p>A very tiny number. It is used for approximating zero.</p>
	 * <!-- end-UML-doc -->
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	private final double tiny = 1.0E-30;

	/** 
	 * <!-- begin-UML-doc -->
	 * <p>The bounds of integration. It is denoted "regn" in the Webnote. It is a vector that describes the lower-left and upper-right coordinates of the Cartesian box that defines the volume of the region. The first n values are the lower-left coordinates and the last n values are the upper-right coordinates, creating a vector of total length 2n where n is the dimensionality of the space.</p>
	 * <!-- end-UML-doc -->
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	private double[] bounds = null;
	/** 
	 * <!-- begin-UML-doc -->
	 * <p>The number of dimensions in the current volume of integration.</p>
	 * <!-- end-UML-doc -->
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	private int numberOfDimensions;
	// All of the other variables from the listing. They are undocumented
	// because I don't know what they do.
	private int i, it, j, k, mds, nd, ndo, ng, npg;
	private double calls, dv2g, dxg, f, f2, f2b, fb, rc, ti;
	private double tsi, wgt, xjac, xn, xnd, xo, schi, si, swgt;
	private int[] ia, kg;
	private double[] dt, dx, r, x, xin;
	// d, di and xi were originally two dimensional matrices, but I have reduced
	// them to 1D for this implementation.
	private double[] d, di, xi;

	/** 
	 * <!-- begin-UML-doc -->
	 * <p>The constructor. It initialized the random number generator and some sets some default values.</p>
	 * <!-- end-UML-doc -->
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public VEGASMonteCarloIntegrator() {
		// begin-user-code

		// Allocate the arrays
		ia = new int[maxDimensions];
		kg = new int[maxDimensions];
		dt = new double[maxDimensions];
		dx = new double[maxDimensions];
		r = new double[maxAxialDivisions];
		x = new double[maxDimensions];
		xin = new double[maxAxialDivisions];
		// Allocate the grids
		d = new double[maxAxialDivisions * maxDimensions];
		di = new double[maxAxialDivisions * maxDimensions];
		xi = new double[maxDimensions * maxAxialDivisions];

		// Setup the random number generator
		randomNumber = new CMWCRNG((int) System.currentTimeMillis());
		//RNG for tests - randomNumber = new CMWCRNG(12345678);

		// end-user-code
	}

	/** 
	 * <!-- begin-UML-doc -->
	 * <p>This operation sets the mode in which the integrator should run. The VEGAS algorithm typically runs in one of four modes, described in this implementation by the NEW_GRID, NEW_ANSWERS, REUSE_ANSWERS and CONTINUE attributes of this class. Unlike the original implementation of VEGAS, this version will only work for those four mode values and the integrate() operation will return 0.0 if the mode is any other value. The definitions of each mode are available in the respective mode's documentation.</p><p>This mode should be called immediately before integrate() and after the region, number of calls and number of iterations have been set.</p>
	 * <!-- end-UML-doc -->
	 * @param mode
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public void setMode(int mode) {
		// begin-user-code

		// Make sure the mode is valid
		if (mode >= 0 && mode <= 3) {
			init_mode = mode;

			numberOfDimensions = bounds.length / 2;
			// Normal entry. Enter here on a cold start.
			if (init_mode <= 0) {
				mds = ndo = 1;
				// Change to mds=0 to disable stratiﬁed sampling, i.e., use
				// importance sampling only.
				for (j = 0; j < numberOfDimensions; j++)
					xi[j * maxAxialDivisions] = 1.0;
			}
			if (init_mode <= 1)
				si = swgt = schi = 0.0;
			// Enter here to inherit the grid from a previous call, but not its
			// answers.
			if (init_mode <= 2) {
				// Enter here to inherit the previous grid and its
				// answers.
				nd = maxAxialDivisions;
				ng = 1;
				// Set up for stratiﬁcation.
				if (mds != 0) {
					ng = (int) Math.pow(maxCalls / 2.0 + 0.25,
							1.0 / numberOfDimensions);
					mds = 1;
					if ((2 * ng - maxAxialDivisions) >= 0) {
						mds = -1;
						npg = ng / maxAxialDivisions + 1;
						nd = ng / npg;
						ng = npg * nd;
					}
				}
				for (k = 1, i = 0; i < numberOfDimensions; i++)
					k *= ng;
				npg = Math.max(maxCalls / k, 2);
				calls = ((double) npg) * ((double) k);

				dxg = 1.0 / ng;
				for (dv2g = 1, i = 0; i < numberOfDimensions; i++)
					dv2g *= dxg;
				dv2g = Math.sqrt(calls * dv2g) / npg / npg / (npg - 1.0);
				xnd = nd;
				dxg *= xnd;
				xjac = 1.0 / calls;
				for (j = 0; j < numberOfDimensions; j++) {
					dx[j] = bounds[j + numberOfDimensions] - bounds[j];
					xjac *= dx[j];
				}
				if (nd != ndo) {
					// Do binning if necessary.
					for (i = 0; i < Math.max(nd, ndo); i++)
						r[i] = 1.0;
					for (j = 0; j < numberOfDimensions; j++)
						rebin(ndo / xnd, j);
					ndo = nd;
				}

				// Print diagnostic information
				// if (nprn >= 0) {
				System.out.println("Input parameters for VEGAS");
				System.out.println("\tndim= " + numberOfDimensions);
				System.out.println("\tncall= " + calls);
				System.out.println("\tit= " + it);
				System.out.println("\titmx= " + maxIterations);
				System.out.println("\tALPH=" + alpha);
				System.out.println("\tmds=" + mds);
				System.out.println("\tnd=" + nd);
				for (j = 0; j < numberOfDimensions; j++) {
					System.out.println("\tx1[" + j + "]= " + bounds[j] + " xu["
							+ j + "]= " + bounds[j + numberOfDimensions]);
				}
				// }
			}

		}

		// end-user-code
	}

	/** 
	 * <!-- begin-UML-doc -->
	 * <p>This operation integrates a function over the integration region. It should be called after the other operations are called to configure the integrator; setMaxIterations(), setIntegrationRegion(), setMode(), setup(), etc.</p>
	 * <!-- end-UML-doc -->
	 * @param function <p>The function that should be integrated over the volume. This function must implemented the IFunction interface. The IFunction.evaluate(x,weight) operation is called by the integrator.</p>
	 * @return <p>An array of three doubles where the first is the value of the integral, the second is the standard deviation and the third is the value of "chi squared per degree of freedom." This final value can be used to determine the quality of the of the integration: if it is much larger than one, then the value of the integration is bad.</p>
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public double[] integrate(IFunction function) {
		// begin-user-code

		// Local Declarations
		int index = 0, secondIndex = 0;
		double integral = 0.0, standardDeviation = 0.0, chiSquaredA = 0.0;
		double[] retArray = new double[3];

		for (it = 0; it < maxIterations; it++) {
			// Main iteration loop. Can enter here (init
			// 3) to do an additional itmx iterations with
			// all other parameters unchanged.
			ti = tsi = 0.0;
			for (j = 0; j < numberOfDimensions; j++) {
				kg[j] = 1;
				for (i = 0; i < nd; i++) {
					index = i * maxDimensions + j; // DONE!
					d[index] = di[index] = 0.0;
				}
			}
			for (;;) {
				fb = f2b = 0.0;
				for (k = 0; k < npg; k++) {
					wgt = xjac;
					for (j = 0; j < numberOfDimensions; j++) {
						xn = (kg[j] - randomNumber.getNextDouble()) * dxg + 1.0;
						ia[j] = Math.max(Math.min((int) xn, maxAxialDivisions),
								1);
						if (ia[j] > 1) {
							index = j * maxAxialDivisions + (ia[j] - 1);
							secondIndex = j * maxAxialDivisions + (ia[j] - 2);
							xo = xi[index] - xi[secondIndex];
							rc = xi[secondIndex] + (xn - ia[j]) * xo;
						} else {
							index = j * maxAxialDivisions + (ia[j] - 1);
							xo = xi[index];
							rc = (xn - ia[j]) * xo;
						}
						x[j] = bounds[j] + rc * dx[j];
						wgt *= xo * xnd;
					}
					f = wgt * function.evaluate(x, wgt);
					f2 = f * f;
					fb += f;
					f2b += f2;
					for (j = 0; j < numberOfDimensions; j++) {
						index = (ia[j] - 1) * maxDimensions + j;
						di[index] += f;
						if (mds >= 0)
							d[index] += f2;
					}
				}
				f2b = Math.sqrt(f2b * npg);
				f2b = (f2b - fb) * (f2b + fb);
				if (f2b <= 0.0)
					f2b = tiny;
				ti += fb;
				tsi += f2b;
				if (mds < 0) {
					// Use stratiﬁed sampling.
					for (j = 0; j < numberOfDimensions; j++) {
						index = (ia[j] - 1) * maxDimensions + j;
						d[index] += f2b;
					}
				}
				for (k = numberOfDimensions - 1; k >= 0; k--) {
					kg[k] %= ng;
					if (++kg[k] != 1)
						break;
				}
				if (k < 0)
					break;
			}
			tsi *= dv2g;
			// Compute final results for this iteration.
			wgt = 1.0 / tsi;
			si += wgt * ti;
			schi += wgt * ti * ti;
			swgt += wgt;
			integral = si / swgt;
			chiSquaredA = (schi - si * integral) / (it + 0.0001);
			if (chiSquaredA < 0.0)
				chiSquaredA = 0.0;
			standardDeviation = Math.sqrt(1.0 / swgt);
			tsi = Math.sqrt(tsi);

			// Print some debug information
			System.out.println("Iteration # " + (it + 1) + ":");
			System.out.println("\tIntegral = " + ti + " +/- " + tsi);

			// It is not clear to me how to fix this comment block because it
			// prints data for each axis based on the value of nprn, which I
			// have eliminated. I think I just need to adjust the loop bounds
			// and print everything.
			// if (nprn != 0) {
			// for (j = 0; j < numberOfDimensions; j++) {
			// System.out.println(" DATA FOR axis " + j);
			// System.out.println("X		delta i		X		delta i");
			// System.out.println("X		deltai");
			// // for (i=nprn/2;i<nd-2;i += nprn+2) {
			// // System.out.println(setw(8) << xi[j][i] << setw(12) <<
			// // di[i][j];
			// // System.out.println(setw(12) << xi[j][i+1] << setw(12) <<
			// // di[i+1][j];
			// // System.out.println(setw(12) << xi[j][i+2] << setw(12) <<
			// // di[i+2][j];
			// // System.out.println(endl;
			// // }
			// }
			// }

			// Reﬁne the grid. Consult references to understand the subtlety of
			// this
			// procedure. The reﬁne- ment is damped, to avoid rapid, destabiliz-
			// ing
			// changes, and also compressed in range by the exponent ALPH.
			for (j = 0; j < numberOfDimensions; j++) {
				xo = d[j];
				xn = d[1 * maxDimensions + j];
				d[j] = (xo + xn) / 2.0;
				dt[j] = d[j];
				for (i = 2; i < nd; i++) {
					rc = xo + xn;
					xo = xn;
					xn = d[i * maxDimensions + j];
					index = (i - 1) * maxDimensions + j;
					d[index] = (rc + xn) / 3.0;
					dt[j] += d[index];
				}
				index = (nd - 1) * maxDimensions + j;
				d[index] = (xo + xn) / 2.0;
				dt[j] += d[index];
			}

			for (j = 0; j < numberOfDimensions; j++) {
				rc = 0.0;
				for (i = 0; i < nd; i++) {
					index = i * maxDimensions + j;
					if (d[index] < tiny)
						d[index] = tiny;
					r[i] = Math.pow((1.0 - d[index] / dt[j])
							/ (Math.log(dt[j]) - Math.log(d[index])), alpha);
					rc += r[i];
				}
				rebin(rc / xnd, j);
			}
		}

		// Setup the array of return values
		retArray[0] = integral;
		retArray[1] = standardDeviation;
		retArray[2] = chiSquaredA;

		System.out.println("All iterations: ");
		System.out
				.println("\tIntegral =" + integral + "+-" + standardDeviation);
		System.out.println("\tchi**2/IT n =" + chiSquaredA);

		return retArray;

		// end-user-code
	}

	/** 
	 * <!-- begin-UML-doc -->
	 * <p>This operation rebins a vector of densities contained in the j-th row of xi into new bins defined by a vector r.</p>
	 * <!-- end-UML-doc -->
	 * @param rc
	 * @param j
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	private void rebin(double rc, int j) {
		// begin-user-code

		// Local Declarations
		int index = 0;

		int i, k = 0;
		double dr = 0.0, xn = 0.0, xo = 0.0;
		for (i = 0; i < nd - 1; i++) {
			while (rc > dr)
				dr += r[(++k) - 1];
			index = j * maxAxialDivisions + (k - 2);
			if (k > 1)
				xo = xi[index];
			index = j * maxAxialDivisions + (k - 1);
			xn = xi[index];
			dr -= rc;
			xin[i] = xn - (xn - xo) * dr / r[k - 1];
		}
		for (i = 0; i < nd - 1; i++)
			xi[j * maxAxialDivisions + i] = xin[i];
		index = j * maxAxialDivisions + (nd - 1);
		xi[index] = 1.0;

		// end-user-code
	}

	/** 
	 * <!-- begin-UML-doc -->
	 * <p>This operation sets the maximum number of times the function will be called during a given iteration.</p>
	 * <!-- end-UML-doc -->
	 * @param maxNumberOfCalls <p>The maximum number of calls to the function.</p>
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public void setMaxNumberOfCalls(int maxNumberOfCalls) {
		// begin-user-code

		maxCalls = maxNumberOfCalls;

		// end-user-code
	}

	/** 
	 * <!-- begin-UML-doc -->
	 * <p>This operation sets the maximum number of iterations that will be performed to refine the solution grid during the integration.</p>
	 * <!-- end-UML-doc -->
	 * @param maxNumberOfIterations <p>The maximum number of iterations.</p>
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public void setMaxNumberOfIterations(int maxNumberOfIterations) {
		// begin-user-code

		maxIterations = maxNumberOfIterations;

		// end-user-code
	}

	/** 
	 * <!-- begin-UML-doc -->
	 * <p>This operation sets the bounds of the integration.</p>
	 * <!-- end-UML-doc -->
	 * @param volume <p>A vector that describes the lower-left and upper-right coordinates of the Cartesian box that defines the volume of the region. The first n values are the lower-left coordinates and the last n values are the upper-right coordinates, creating a vector of total length 2n where n is the dimensionality of the space.</p>
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public void setIntegrationRegion(double[] volume) {
		// begin-user-code

		// Set the bounds of the integration
		bounds = volume;

		// end-user-code
	}
}