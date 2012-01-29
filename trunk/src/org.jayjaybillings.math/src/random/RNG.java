/**
 * 
 */
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