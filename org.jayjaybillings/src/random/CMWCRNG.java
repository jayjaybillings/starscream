/**
 * 
 */
package random;

import java.util.Random;

/** 
 * <!-- begin-UML-doc -->
 * <p>A complimentary multiply-with-carry random number generator with a period of on the order of 2^131104. Based on the work from George Marsaglia (http://school.anhb.uwa.edu.au/personalpages/kwessen/shared/Marsaglia03.html), the implementation provided on Wikipedia, (http://en.wikipedia.org/wiki/Multiply-with-carry), and with the weight "a" taken from Numerical Recipes, 3rd. Edition. This random number generator has been verified to pass the Diehard tests and other tests in the Dieharder library, (http://www.phy.duke.edu/~rgb/General/dieharder.php).</p>
 * <!-- end-UML-doc -->
 * @author jaybilly
 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
 */
public class CMWCRNG extends RNG {

	/** 
	 * <!-- begin-UML-doc -->
	 * <p>State vector index for the generator</p>
	 * <!-- end-UML-doc -->
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	private int index;
	/** 
	 * <!-- begin-UML-doc -->
	 * <p>The addend used for the generator.</p>
	 * <!-- end-UML-doc -->
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	private long addend;

	/** 
	 * <!-- begin-UML-doc -->
	 * <p>The constructor</p>
	 * <!-- end-UML-doc -->
	 * @param seed
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public CMWCRNG(int seed) {
		// begin-user-code

		// Call the super constructor
		super(seed);

		// Create an instance of java.util.Random for initializing the state
		// vector
		Random stdRng = null;
		// Use the user suplied seed
		if (seed != 0)
			stdRng = new Random(seed);
		else
			// Otherwise use the clock value
			stdRng = new Random(System.nanoTime());

		// Set the state size
		stateSize = 4096;

		// Initialize the state vector
		stateVector = new long[stateSize];
		for (int i = 0; i < stateSize; i++) {
			stateVector[i] = stdRng.nextLong();
		}

		return;
		// end-user-code
	}

	/** 
	 * <!-- begin-UML-doc -->
	 * <p>The copy constructor</p>
	 * <!-- end-UML-doc -->
	 * @param original
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public CMWCRNG(RNG original) {
		// begin-user-code

		// Call the super constructor
		super(original);

		// Copy the state size and values
		original.stateSize = stateSize;
		for (int i = 0; i < stateSize; i++) {
			original.stateVector[i] = stateVector[i];
		}

		// end-user-code
	}

	/** 
	 * (non-Javadoc)
	 * @see RNG#getNextInt()
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public int getNextInt() {
		// begin-user-code
		return (int) getNextLong();
		// end-user-code
	}

	/** 
	 * <!-- begin-UML-doc -->
	 * <p>This operation returns a randomly generated 64-bit "long" integer.</p>
	 * <!-- end-UML-doc -->
	 * @return <p>The long integer.</p>
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	protected long getNextLong() {
		// begin-user-code
		// Local Declarations
		long t, multiplier = 4294957665L;
		index = stateSize - 1;
		long nextLong, base = 0xfffffffe;

		// Generate the random number
		index = (index + 1) & (stateSize - 1);
		t = multiplier * stateVector[index] + addend;
		addend = (t >>> 32);
		nextLong = t + addend;
		if (nextLong < addend) {
			nextLong++;
			addend++;
		}

		// Return the number after storing it
		return (stateVector[index] = base - nextLong);

		// end-user-code
	}

	/** 
	 * (non-Javadoc)
	 * @see RNG#getNextDouble()
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public double getNextDouble() {
		// begin-user-code
		return 5.42101086242752217E-20 * getNextLong();
		// end-user-code
	}

	/** 
	 * (non-Javadoc)
	 * @see RNG#getNextFloat()
	 * @generated "UML to Java (com.ibm.xtools.transform.uml2.java5.internal.UML2JavaTransform)"
	 */
	public float getNextFloat() {
		// begin-user-code
		return (float) getNextDouble();
		// end-user-code
	}
}