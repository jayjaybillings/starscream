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
		long startTime = 0, quittingTime = 0,
		     runTime = 3600000; // 10 hours in milliseconds
		int value = 0;
		
		// Setup starting time
		startTime = currentDate.getTime();
		quittingTime = startTime + runTime;
				
		// Loop until the max time is met
		while(currentDate.getTime() < quittingTime) {
			// Get the next random integer value
			value = cmwcGenerator.getNextInt();
			// Write the value to the output stream in binary
			System.out.write(value);
		}

		return;
		// end-user-code
	}
}