module str

    !this module has some string manipulation commands 
    public
    contains

    pure function tok_count(text)
        !counts number of tokens in a string
        character(len =  *), intent(in) ::  text
        integer :: tok_count
        integer :: i, j
        logical :: in_tok

        j = len(trim(adjustl(text)))
        in_tok = .false.
        tok_count = 0
        do i = 1, j
            !This checks if it is a white space character which is the delimiter
            if(trim(adjustl(text(i:i))) == ' ') then 
                !If previously we were in token and the current character is the delimiter 
                !Then we are no longer in the token
                if(in_tok) in_tok = .false.

            !If the character isn't a white space character and we previously weren't in the token then set in_tok 
            !to true and increment token count
            else if(.not.in_tok) then 
                in_tok = .true.
                tok_count = tok_count + 1
            end if
        end do
        return
    end function tok_count

       subroutine to_lower(str)
     character(*), intent(in out) :: str
     integer :: i
 
     do i = 1, len(str)
       select case(str(i:i))
         case("A":"Z")
           str(i:i) = achar(iachar(str(i:i))+32)
       end select
     end do  
   end subroutine to_lower

end module str
