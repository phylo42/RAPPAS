/*
 * Copyright (C) 2018 yann
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package etc.exceptions;

import core.DNAStates;
import core.DNAStatesShifted;
import core.States;

/**
 * exception to report characters not compatible with IUPAC code
 * @author yann
 */
public class NonIUPACStateException extends Exception {
    
    char c;

    public NonIUPACStateException() {super();}
    
    public NonIUPACStateException(States s, char c) {
        super("The non-IUPAC state "+c+" is not allowed.");
        this.c=c;
        if ( (s instanceof DNAStates) || (s instanceof DNAStatesShifted)) {
            System.err.println("You selected a nucleotide analysis (-s 'nucl'), but this state is not in IUPAC allowed nucleotides.");
        } else {
            System.err.println("You selected an amino acid analysis (-s 'amino'), but this state is not in IUPAC allowed amino acids.");
        }
    }
    
    public char getFaultyChar() {
        return c;
    }
}
