/**
 * 
 */
package integration;

import common.IFunction;

/** 
 * <!-- begin-UML-doc -->
 * <p>A 3D sinc IFunction used in the test. f(x,y,z) = sin(x)/x + sin(y)/y + sin(z)/z over the region -20.0 &lt;= {x,y,z} &lt;= 20.0. The value of the integral is 14863.1 according to Wolfram Alpha.</p>
 * <!-- end-UML-doc -->
 * @author jaybilly
 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
 */
public class ThreeDSinc implements IFunction {
	/** 
	 * (non-Javadoc)
	 * @see IFunction#evaluate(double[] x)
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public double evaluate(double[] x) {
		// begin-user-code

		// Local Declarations
		double x_i = x[0], x_j = x[1], x_k = x[2];

		return (Math.sin(x_i) / x_i) + (Math.sin(x_j) / x_j)
				+ (Math.sin(x_k) / x_k);
		// end-user-code
	}

	/** 
	 * (non-Javadoc)
	 * @see IFunction#evaluate(double[] x, double weight)
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public double evaluate(double[] x, double weight) {
		// begin-user-code
		return evaluate(x);
		// end-user-code
	}

	/** 
	 * (non-Javadoc)
	 * @see IFunction#getDimensionality()
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public int getDimensionality() {
		// begin-user-code
		return 3;
		// end-user-code
	}
}