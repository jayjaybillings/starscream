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
package common;

/** 
 * <!-- begin-UML-doc -->
 * <p>The IFunction interface represents a mathematical function. </p>
 * <!-- end-UML-doc -->
 * @author jaybilly
 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
 */
public interface IFunction {
	/** 
	 * <!-- begin-UML-doc -->
	 * <p>This operation evaluates the function for the vector x of length n. If n is greater than the dimensionality d of the function, then the function is evaluated using the first d values of x. However, if n is less than d, then function will return 0.0.</p>
	 * <!-- end-UML-doc -->
	 * @param x <p>A vector of length n = d where d is the function's dimensionality at which the function is evaluated.</p>
	 * @return <p>The value of the function at x.</p>
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public double evaluate(double[] x);

	/** 
	 * <!-- begin-UML-doc -->
	 * <p>This operation works the same as evaluate(x), but with an additional weight argument than can be used to integrate another function g(x) that resembles f.</p>
	 * <!-- end-UML-doc -->
	 * @param x <p>A vector of length n = d where d is the function's dimensionality at which the function is evaluated.</p>
	 * @param weight <p>A weight w that can be used to integrate a function g(x) that resembles this function using the relation integral(g) = sum[0,i](w_i*g(x)).</p>
	 * @return <p>The value of the function at x.</p>
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public double evaluate(double[] x, double weight);

	/** 
	 * <!-- begin-UML-doc -->
	 * <p>This operation returns the dimensionality of the function.</p>
	 * <!-- end-UML-doc -->
	 * @return
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public int getDimensionality();
}