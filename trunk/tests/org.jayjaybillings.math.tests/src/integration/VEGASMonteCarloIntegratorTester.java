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
/**
 * 
 */
package integration;

import static org.junit.Assert.*;
import org.junit.Test;

import common.IFunction;

/** 
 * <!-- begin-UML-doc -->
 * <p>This class is responsible for testing the VEGASMonteCarloIntegrator. It integrates a set of functions, most of which have known analytic integrals, and compares the answer from the integrator to a gold standard, either the analytical solution or the best obtainable numerical solution.</p>
 * <!-- end-UML-doc -->
 * @author jaybilly
 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
 */
public class VEGASMonteCarloIntegratorTester {
	/** 
	 * <!-- begin-UML-doc -->
	 * <p>An instance of the VEGASMonteCarloIntegrator that is used for the tests. The same integrator is used for all of the tests.</p>
	 * <!-- end-UML-doc -->
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	private VEGASMonteCarloIntegrator integrator;

	/** 
	 * <!-- begin-UML-doc -->
	 * <p>This operation integrates a simple one dimensional quadratic function f(x) = ax^2+b^x+c over the region -100.0 &lt; x &lt; 100.0. For a = 24.3, b = 81.0 and c = 5.0, the value of the integral is 1.6201*10^7.</p>
	 * <!-- end-UML-doc -->
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	@Test
	public void check1DQuadratic() {
		// begin-user-code

		// Local Declarations
		IFunction function = new OneDQuadratic();
		double[] volume = new double[2], answer = new double[3];
		double analyticAnswer = 1.6201E7;

		// Setup the integrator. Manually testing this function found that it
		// gets the correct answer with a low value of (Chi^2)*a, normally
		// around 2.8, for the parameters below. Since this is a very simple and
		// smooth function without significant features, the number of calls to
		// the function should be more significant than the number of
		// iterations. The integration should be accurate up to 3%.
		integrator = new VEGASMonteCarloIntegrator();
		volume[0] = -100.0;
		volume[1] = 100.0;
		integrator.setIntegrationRegion(volume);
		integrator.setMaxNumberOfCalls(15000);
		integrator.setMaxNumberOfIterations(2);
		integrator.setMode(0);

		// Integrate the function
		answer = integrator.integrate(function);

		// Check the value of the integral. It should be accurate up to 3% for
		// this function.
		assertEquals(analyticAnswer, answer[0], analyticAnswer / 33.0);

		// Write the output for diagnostics
		printIntegralValue(answer[0], answer[1], answer[2], "OneDQuadratic");

		return;

		// end-user-code
	}

	/** 
	 * <!-- begin-UML-doc -->
	 * <p>This operation integrates the UnitSphere function to calculate its volume in the region defined by the unit cube (-1&lt;={x,y,z}&lt;=1). The volume of a unit sphere (R=1) is (4/3)*pi.</p>
	 * <!-- end-UML-doc -->
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	@Test
	public void checkUnitSphereVolume() {
		// begin-user-code

		// Local Declarations
		IFunction function = new UnitSphere();
		double[] volume = new double[6], answer = new double[3];
		double analyticAnswer = (4.0 / 3.0) * Math.PI; // 4.19

		// Setup the integrator. The integration should be accurate up to 3%.
		integrator = new VEGASMonteCarloIntegrator();
		volume[0] = -1.0;
		volume[1] = -1.0;
		volume[2] = -1.0;
		volume[3] = 1.0;
		volume[4] = 1.0;
		volume[5] = 1.0;
		integrator.setIntegrationRegion(volume);
		integrator.setMaxNumberOfCalls(10000);
		integrator.setMaxNumberOfIterations(2);
		integrator.setMode(0);

		// Integrate the function
		answer = integrator.integrate(function);

		// Check the value of the integral. It should be accurate up to 3% for
		// this function.
		assertEquals(analyticAnswer, answer[0], analyticAnswer / 33.0);

		// Write the output for diagnostics
		printIntegralValue(answer[0], answer[1], answer[2], "Unit Sphere");

		// end-user-code
	}

	/** 
	 * <!-- begin-UML-doc -->
	 * <p>This operation integrates a 3D sinc, f(x,y,z) = sin(x)/x + sin(y)/y + sin(z)/z, over the region -20.0 &lt;= {x,y,z} &lt;= 20.0. The value of the integral is 14863.1 according to Wolfram Alpha.</p>
	 * <!-- end-UML-doc -->
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	@Test
	public void check3DSinc() {
		// begin-user-code

		// Local Declarations
		IFunction function = new ThreeDSinc();
		double[] volume = new double[6], answer = new double[3];
		double analyticAnswer = 14863.1; // According to Wolfram Alpha

		// Setup the integrator. The integration should be accurate up to 3%.
		integrator = new VEGASMonteCarloIntegrator();
		volume[0] = -20.0;
		volume[1] = -20.0;
		volume[2] = -20.0;
		volume[3] = 20.0;
		volume[4] = 20.0;
		volume[5] = 20.0;
		integrator.setIntegrationRegion(volume);
		integrator.setMaxNumberOfCalls(10000);
		integrator.setMaxNumberOfIterations(3);
		integrator.setMode(0);

		// Integrate the function
		answer = integrator.integrate(function);

		// Check the value of the integral. It should be accurate up to 3% for
		// this function.
		assertEquals(analyticAnswer, answer[0], analyticAnswer / 33.0);

		// Write the output for diagnostics
		printIntegralValue(answer[0], answer[1], answer[2], "3D Sinc Function");

		// end-user-code
	}

	/** 
	 * <!-- begin-UML-doc -->
	 * <p>This operation prints the value of an integration to the screen. It is simply a utility function to avoid repeating the dump to stdout.</p>
	 * <!-- end-UML-doc -->
	 * @param value <p>The value of the integral.</p>
	 * @param stdDev <p>The standard deviation</p>
	 * @param chi2a <p>(Chi^2)*a</p>
	 * @param testName <p>The name of the test (OneDQuadratic, etc.).</p>
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	private void printIntegralValue(double value, double stdDev, double chi2a,
			String testName) {
		// begin-user-code

		// Dump to stdout
		System.out.println("Integration results for " + testName + ":");
		System.out.println("\tIntegral value = " + value);
		System.out.println("\tStandard Deviation = " + stdDev);
		System.out.println("\t(Chi^2)*a = " + chi2a + "\n");

		return;
		// end-user-code
	}
}