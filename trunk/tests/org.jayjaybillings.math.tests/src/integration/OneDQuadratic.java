/**
 * 
 */
package integration;

import common.IFunction;

/** 
 * <!-- begin-UML-doc -->
 * <p>A 1D quadratic IFunction used in the test. f(x) = ax^2+b^x+c over the region -100.0 &lt; x &lt; 100.0. For a = 24.3, b = 81.0 and c = 5.0, the value of the integral is 1.6201*10^7.</p>
 * <!-- end-UML-doc -->
 * @author jaybilly
 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
 */
public class OneDQuadratic implements IFunction {
	/** 
	 * (non-Javadoc)
	 * @see IFunction#evaluate(double[] x)
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public double evaluate(double[] x) {
		// begin-user-code
		return 24.3 * x[0] * x[0] + 81 * x[0] + 5.0;
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
		return 1;
		// end-user-code
	}
}