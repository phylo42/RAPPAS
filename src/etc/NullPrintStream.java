/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package etc;

/**
 *
 * @author ben
 */
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;

public class NullPrintStream extends PrintStream {

  public NullPrintStream() {
    super(new NullByteArrayOutputStream());
  }

  private static class NullByteArrayOutputStream extends ByteArrayOutputStream {

    @Override
    public void write(int b) {
      // do nothing
    }

    @Override
    public void write(byte[] b, int off, int len) {
      // do nothing
    }

    @Override
    public void writeTo(OutputStream out) throws IOException {
      // do nothing
    }

  }

}
