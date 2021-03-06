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
package random;

import java.util.Date;

/** 
 * <!-- begin-UML-doc -->
 * <p>This class generates output for the Dieharder RNG testing suite using the CMWC-RNG.</p>
 * <!-- end-UML-doc -->
 * @author jaybilly
 */
public class CMWCRNGOutputGenerator {
	/** 
	 * <!-- begin-UML-doc -->
	 * <p>The&nbsp;all-powerful&nbsp;main()!&nbsp;This&nbsp;operation&nbsp;generates&nbsp;the&nbsp;output&nbsp;and&nbsp;writes&nbsp;it&nbsp;to&nbsp;a&nbsp;binary&nbsp;stream&nbsp;on&nbsp;stdout.</p>
	 * <!-- end-UML-doc -->
	 * @param args
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public static void main(String[] args) {
		// begin-user-code

		// Local Declarations - all times in milliseconds
		CMWCRNG cmwcGenerator = new CMWCRNG((int) System.nanoTime());
		Date currentDate = new Date();
		long startTime = 0, quittingTime = 0, runTime = 3600000; // 10 hours in milliseconds
		int value = 0;

		// Setup starting time
		startTime = currentDate.getTime();
		quittingTime = startTime + runTime;

		// Loop until the max time is met
		while (currentDate.getTime() < quittingTime) {
			// Get the next random integer value
			value = cmwcGenerator.getNextInt();
			// Write the value to the output stream in binary
			System.out.write(value);
		}

		return;
		// end-user-code
	}
}