/**
 * 
 */
package integration;

import common.IFunction;

/** 
 * <!-- begin-UML-doc -->
 * <p>An IFunction for the unit sphere (R=1) used for testing the VEGASMonteCarloIntegrator. The volume of a unit sphere is (4/3)*pi.</p>
 * <!-- end-UML-doc -->
 * @author jaybilly
 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
 */
public class UnitSphere implements IFunction {
	/** 
	 * (non-Javadoc)
	 * @see IFunction#evaluate(double[] x)
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public double evaluate(double[] x) {
		// begin-user-code

		// Local Declarations
		double radius = Math.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);

		// Return 1.0 if the radius is less than one, zero otherwise 
		return (radius <= 1.0) ? 1.0 : 0.0;
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