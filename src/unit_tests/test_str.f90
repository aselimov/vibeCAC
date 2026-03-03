program test_str
    use str, only: tok_count, to_lower
    use assertions, only: assert_equal_int, assert_equal_str
    implicit none

    integer :: failures

    failures = 0

    call test_tok_count_empty_string(failures)
    call test_tok_count_whitespace_only(failures)
    call test_tok_count_trimmed_multi_tokens(failures)
    call test_tok_count_many_tokens(failures)
    call test_to_lower_mixed_case(failures)
    call test_to_lower_keeps_non_letters(failures)

    if (failures > 0) then
        print '(A,1X,I0)', 'Unit tests failed:', failures
        stop 1
    end if

    print '(A)', 'All str unit tests passed.'

contains

    subroutine test_tok_count_empty_string(failures)
        integer, intent(inout) :: failures

        ! Arrange / Act / Assert
        call assert_equal_int(tok_count(''), 0, 'empty string should have zero tokens', failures)
    end subroutine test_tok_count_empty_string

    subroutine test_tok_count_whitespace_only(failures)
        integer, intent(inout) :: failures

        ! Arrange / Act / Assert
        call assert_equal_int(tok_count('      '), 0, 'whitespace-only input should have zero tokens', failures)
    end subroutine test_tok_count_whitespace_only

    subroutine test_tok_count_trimmed_multi_tokens(failures)
        integer, intent(inout) :: failures

        ! Arrange / Act / Assert
        call assert_equal_int(tok_count('  run   2500   '), 2, 'leading/trailing and repeated spaces should be ignored', failures)
    end subroutine test_tok_count_trimmed_multi_tokens

    subroutine test_tok_count_many_tokens(failures)
        integer, intent(inout) :: failures

        ! Arrange / Act / Assert
        call assert_equal_int(tok_count('a b c d e f g h i j'), 10, 'many valid tokens should be counted correctly', failures)
    end subroutine test_tok_count_many_tokens

    subroutine test_to_lower_mixed_case(failures)
        integer, intent(inout) :: failures
        character(len=16) :: text

        ! Arrange
        text = 'ReAd_DaTa'

        ! Act
        call to_lower(text)

        ! Assert
        call assert_equal_str(trim(text), 'read_data', 'to_lower should convert uppercase letters', failures)
    end subroutine test_to_lower_mixed_case

    subroutine test_to_lower_keeps_non_letters(failures)
        integer, intent(inout) :: failures
        character(len=16) :: text

        ! Arrange
        text = 'TEMP_300-STEP2'

        ! Act
        call to_lower(text)

        ! Assert
        call assert_equal_str(trim(text), 'temp_300-step2', 'to_lower should leave digits/symbols unchanged', failures)
    end subroutine test_to_lower_keeps_non_letters

end program test_str
