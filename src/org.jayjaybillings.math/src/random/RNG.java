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

package random;

/** 
 * <!-- begin-UML-doc -->
 * <p>The RNG class is the abstract base class for all random number generators in the org.jayjaybillings.math package. It uses java.util.Random to set its initial state.</p>
 * <!-- end-UML-doc -->
 * @author jaybilly
 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
 */
public abstract class RNG {
	/** 
	 * <!-- begin-UML-doc -->
	 * <p>The state vector of the RNG.</p>
	 * <!-- end-UML-doc -->
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	protected long[] stateVector = null;

	/** 
	 * <!-- begin-UML-doc -->
	 * <p>The size of the state vector.</p>
	 * <!-- end-UML-doc -->
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	protected int stateSize;

	/** 
	 * <!-- begin-UML-doc -->
	 * <p>The constructor</p>
	 * <!-- end-UML-doc -->
	 * @param seed <p>The seed value with which the RNG should be initialized.</p>
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public RNG(int seed) {
		// begin-user-code
		// end-user-code
	}

	/** 
	 * <!-- begin-UML-doc -->
	 * <p>A copy constructor. This constructor should be used in cases where a copy of the RNG that produces the exact same numbers as the original is needed. This is particularly useful for testing and debuging.</p>
	 * <!-- end-UML-doc -->
	 * @param original <p>A second RNG whose state should be copied to create a new RNG.</p>
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public RNG(RNG original) {
		// begin-user-code
		// end-user-code
	}

	/** 
	 * <!-- begin-UML-doc -->
	 * <p>This operation returns a randomly generated 32-bit integer.</p>
	 * <!-- end-UML-doc -->
	 * @return <p>The integer</p>
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public abstract int getNextInt();

	/** 
	 * <!-- begin-UML-doc -->
	 * <p>This operation returns a randomly generated 64-bit double-precision floating point number.</p>
	 * <!-- end-UML-doc -->
	 * @return <p>The double-precision floating-point number.</p>
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public abstract double getNextDouble();

	/** 
	 * <!-- begin-UML-doc -->
	 * <p>This operation returns a randomly generated 32-bit single-precision floating point number.</p>
	 * <!-- end-UML-doc -->
	 * @return <p>The single-precision floating-point number.</p>
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public abstract float getNextFloat();
}