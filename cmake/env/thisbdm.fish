# -----------------------------------------------------------------------------
#
# Copyright (C) The BioDynaMo Project.
# All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
#
# See the LICENSE file distributed with this work for details.
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# -----------------------------------------------------------------------------
# Source this script to set up the BioDynaMo build that this script is part of.
#
# Conveniently an alias like this can be defined in your config.fish:
#   alias thisbdm="source path/to/biodynamo/bin/thisbdm.fish -q"
#
# This script is for fish, see thisbdm.sh for bash like shells.
#
# 'thisbdm.sh' author: Fons Rademakers, 16/5/2019
# fish port author: M. Dorukhan Arslan, 26/7/2020

function _drop_from_var
    if test (count $argv) -ne 2
        echo "_drop_from_var: requires 2 arguments but has " (count $argv)
        return 1
    end

    set -l var $argv[1]
    set -l drop $argv[2]

    if  not set -q $var; or test -z "$drop"
        return 0
    end

    if set -l index (contains -i $drop $$var)
        set -e "$var"["$index"]
    end
end

function _bdm_err
    $_bdm_silent; or echo -e "\e[31m$argv\e[0m"
end

function _bdm_ok
    $_bdm_quiet; or echo -e "\e[32m$argv\e[0m"
end

function _thisbdm_cleanup
  functions -e source_thisbdm
  functions -e _drop_from_var
  functions -e _bdm_err
  functions -e _bdm_ok
  set -e _bdm_silent
  set -e _bdm_quiet
  functions -e _thisbdm_cleanup
end

function source_thisbdm
    set -g _bdm_quiet false
    set -g _bdm_silent false
    for option in $argv
        switch "$option"
            case -q --quiet
                set _bdm_quiet true
            case -Q --silent
                set _bdm_quiet true
                set _bdm_silent true
            case \*
                return 1
        end
    end

    set -l old_bdmsys
    if test -n "$BDMSYS"
        set old_bdmsys $BDMSYS
    end

    set -l curr_filename (status --current-filename)
    set -gx BDMSYS (fish -c "cd (dirname $curr_filename)/..; and pwd"); or return 1

    # Clear the env from previously set BioDynaMo paths.
    if test -n "$old_bdmsys"
        _drop_from_var PATH "$old_bdmsys/bin"
        _drop_from_var LD_LIBRARY_PATH "$old_bdmsys/lib"
        _drop_from_var DYLD_LIBRARY_PATH "$old_bdmsys/lib"
        _drop_from_var SHLIB_PATH "$old_bdmsys/lib"
        _drop_from_var MANPATH "$old_bdmsys/man"
        _drop_from_var CMAKE_PREFIX_PATH "$old_bdmsys"
	end

	if test -n "$old_bdmsys"
        _drop_from_var ParaView_DIR "$old_bdmsys/third_party/paraview/lib/cmake/paraview-5.8"
        _drop_from_var ParaView_LIB_DIR "$old_bdmsys/third_party/paraview/lib"
        _drop_from_var PV_PLUGIN_PATH "$old_bdmsys/biodynamo/lib/pv_plugin"
        _drop_from_var PATH "$old_bdmsys/third_party/paraview/bin"
        _drop_from_var Qt5_DIR "$old_bdmsys/third_party/qt/lib/cmake/Qt5"
        _drop_from_var QT_QPA_PLATFORM_PLUGIN_PATH "$old_bdmsys/third_party/qt/plugins"
        _drop_from_var DYLD_LIBRARY_PATH "$old_bdmsys/third_party/paraview/lib"
        _drop_from_var DYLD_LIBRARY_PATH "$old_bdmsys/third_party/qt/lib"
        _drop_from_var LD_LIBRARY_PATH "$old_bdmsys/third_party/paraview/lib"
        _drop_from_var LD_LIBRARY_PATH "$old_bdmsys/third_party/qt/lib"
	end
    #########

    set -l default_manpath
    if test -z "$MANPATH"
        # Grab the default man path before setting the path to avoid duplicates
        if type -qt manpath
            set default_manpath (manpath); or return 1
        else if type -qt man
            set default_manpath (man -w 2> /dev/null); or return 1
        end
    end

    if test -z "$PATH"
        set -gx PATH "$BDMSYS/bin"
    else
        set -gx PATH "$BDMSYS/bin:$PATH"
    end

    if test -z "$LD_LIBRARY_PATH"
        set -gx LD_LIBRARY_PATH "$BDMSYS/lib" # Linux, ELF HP-UX
    else
        set -pgx LD_LIBRARY_PATH "$BDMSYS/lib"
    end

    if test -z "$DYLD_LIBRARY_PATH"
        set -gx DYLD_LIBRARY_PATH "$BDMSYS/lib" # Mac OS X
    else
        set -pgx DYLD_LIBRARY_PATH "$BDMSYS/lib"
    end

    if test -z "$SHLIB_PATH"
        set -gx SHLIB_PATH "$BDMSYS/lib" # legacy HP-UX
    else
        set -pgx SHLIB_PATH "$BDMSYS/lib"
    end

    if test -z "$LIBPATH"
        set -gx LIBPATH "$BDMSYS/lib" # AIX
    else
        set -pgx LIBPATH "$BDMSYS/lib"
    end

    if test -z "$MANPATH"
        set -gx MANPATH "$BDMSYS/man":$default_manpath
    else
        set -pgx MANPATH "$BDMSYS/man"
    end

    ##### Python Specific Configurations #####
    set -gx PYENV_ROOT @pyenvroot@
    if test -z "$PYENV_ROOT"
        set -gx PYENV_ROOT "$HOME/.pyenv"
    end

    set -pgx PATH "$PYENV_ROOT/bin"

    # FIXME Some paths are (ap/pre)pended n times for n calls to thisbdm.*sh
    # due to https://github.com/pyenv/pyenv/issues/969
    pyenv init - | source; or return 1
    pyenv shell @pythonvers@; or return 1

    # Location of jupyter executable (installed with `pip install --user` command)
    if test -n "$PYTHONUSERBASE"
        set -pgx PATH "$PYTHONUSERBASE/.local/bin"
    else
        set -pgx PATH "$HOME/.local/bin"
    end
    set -pgx LD_LIBRARY_PATH "$PYENV_ROOT/versions/@pythonvers@/lib"
    ########

    ##### CMake Specific Configurations #####
    if test -z "$CMAKE_PREFIX_PATH"
        set -gx CMAKE_PREFIX_PATH "$BDMSYS/share/cmake"  # Linux, ELF HP-UX
    else
        set -pgx CMAKE_PREFIX_PATH "$BDMSYS/share/cmake"
    end
    ########

    #### ROOT Specific Configurations ####
    set -l orvers
    set -l crvers

    if test -z "$BDM_ROOT_DIR"; and test -z "$ROOTSYS"
        set -gx BDM_ROOT_DIR "$BDMSYS/third_party/root"
        if not test -d "$BDM_ROOT_DIR"
            _bdm_err "[ERR] We are unable to source ROOT! Please make sure ROOT is installed on your system!"
            _bdm_err "      You can specify manually its location by executing 'set -gx BDM_ROOT_DIR path/to/root'"
            _bdm_err "      before running cmake."
            return 1
        end
    else
        # ROOTSYS has precedence over the BDM_ROOT_DIR custom configuration
        if test -n "$ROOTSYS"
            set orvers "@rootvers@"
            set crvers ("$ROOTSYS"/bin/root-config --version)
            if test "$crvers" = "$orvers"
                set -gx BDM_ROOT_DIR "$ROOTSYS"
            else
                _bdm_err "[ERR] ROOTSYS points to ROOT version $crvers, while BDM was build with version $orvers."
                _bdm_err "      Make sure that ROOTSYS points to the right version of ROOT."
                return 1
            end
        end
    end

    function __bdm_root
        "$BDM_ROOT_DIR"/bin/root -l -e 'cout << "Loading BioDynaMo into ROOT..." << endl; gROOT->LoadMacro("'"$BDMSYS"'/etc/rootlogon.C");' $argv
    end
    funcsave __bdm_root
    . "$BDM_ROOT_DIR/bin/thisroot.fish"
    ########

    #### ParaView Specific Configurations ####
    set -l with_paraview @with_paraview@
    if test "$with_paraview" = 'ON'
        if test -z "$ParaView_DIR"
            set -gx ParaView_DIR $BDMSYS/third_party/paraview
            if not test -d "$ParaView_DIR"
                _bdm_err "[ERR] We are unable to find ParaView! Please make sure it is installed in your system!"
                _bdm_err "      You can specify manually its location by executing 'set -gx ParaView_DIR path/to/paraview'"
                _bdm_err "      together with 'set -gx Qt5_DIR path/to/qt' before running cmake."
                return 1
            end
        end

        if test -z "$ParaView_LIB_DIR"
            set -gx ParaView_LIB_DIR "$ParaView_DIR/lib"
        else
            set -pgx ParaView_LIB_DIR "$ParaView_DIR/lib"
        end

        if test -z "$PV_PLUGIN_PATH"
            set -gx PV_PLUGIN_PATH "$BDMSYS/lib/pv_plugin"
        else
            set -pgx PV_PLUGIN_PATH "$BDMSYS/lib/pv_plugin"
        end

        # We don't add the ParaView site-packages path to PYTHONPATH, because pip in the
        # pyenv environment will not function anymore: ModuleNotFoundError: No module named 'pip._internal'
        alias __bdm_paraview='$ParaView_DIR/bin/paraview'; funcsave __bdm_paraview
        # aliases are just wrapped functions in fish, so they have the desired behavior
        alias __bdm_pvpython='$ParaView_DIR/bin/pvpython'; funcsave __bdm_pvpython
        alias __bdm_pvbatch='$ParaView_DIR/bin/pvbatch'; funcsave __bdm_pvbatch

        if test -z "$LD_LIBRARY_PATH"
            set -gx LD_LIBRARY_PATH "$ParaView_LIB_DIR"
        else
            set -pgx LD_LIBRARY_PATH "$ParaView_LIB_DIR"
        end

        if test -z "$DYLD_LIBRARY_PATH"
            set -gx DYLD_LIBRARY_PATH "$ParaView_LIB_DIR"
        else
            set -pgx DYLD_LIBRARY_PATH "$ParaView_LIB_DIR"
        end
        ########

        #### Qt5 Specific Configurations ####
        if test -z "$Qt5_DIR"
            set -gx Qt5_DIR $BDMSYS/third_party/qt
            if not test -d "$Qt5_DIR"
                _bdm_err "[ERR] We are unable to find Qt! Please make sure it is installed in your system!"
                _bdm_err "      You can specify manually its location by executing 'set -gx Qt5_DIR path/to/qt'"
                _bdm_err "      together with 'set -gx ParaView_DIR path/to/paraview' before running cmake."
                return 1
            end
        end

        if test -z "$QT_QPA_PLATFORM_PLUGIN_PATH"
            set -gx QT_QPA_PLATFORM_PLUGIN_PATH "$Qt5_DIR/plugins"
        else
            set -pgx QT_QPA_PLATFORM_PLUGIN_PATH "$Qt5_DIR/plugins"
        end

        if test -z "$LD_LIBRARY_PATH"
            set -gx LD_LIBRARY_PATH "$Qt5_DIR/lib"
        else
            set -pgx LD_LIBRARY_PATH "$Qt5_DIR/lib"
        end

        if test -z "$DYLD_LIBRARY_PATH"
            set -gx DYLD_LIBRARY_PATH "$Qt5_DIR/lib"
        else
            set -pgx DYLD_LIBRARY_PATH "$Qt5_DIR/lib"
        end
    end
    #######

    # OpenMP
    set -gx OMP_PROC_BIND true

    ###### Platform-specific Configuration
    # Apple specific
    if test (uname) = 'Darwin'
        # Nothing for now
        true
    else # GNU/Linux
        # CentOS specifics
        set -l os_id (grep -oP '(?<=^ID=).+' /etc/os-release | tr -d '"'); or return 1
        if test "$os_id" = 'centos'
            set -gx MESA_GL_VERSION_OVERRIDE "3.3"
            if test -z "$CXX"; and test -z "$CC"
                . scl_source enable devtoolset-7; or return 1
            end

            . /etc/profile.d/modules.sh; or return 1
            module load mpi; or return 1

            # load llvm 6 required for libroadrunner
            if test -d "$BDMSYS"/third_party/libroadrunner
                . scl_source enable llvm-toolset-6.0; or return 1
            end
        end
    end
    #######

    ### Enable commands in child shells (like in bash) ###
    function __bdm_fish_functions
        if test -d "$BDMSYS"
            if test -d "$BDM_ROOT_DIR"
                alias root='__bdm_root'
            end
            if test -d "$ParaView_DIR"
                alias paraview='__bdm_paraview'
                alias pvpython='__bdm_pvpython'
                alias pvbatch='__bdm_pvbatch'
            end
        end
    end
    funcsave __bdm_fish_functions
    set -l marker ' # >>thisbdm<<'
    if test -e $__fish_config_dir/config.fish
            # ensure the above is only called once in config.fish
            sed -i.bak '/^.*'"$marker"'$/,$d' $__fish_config_dir/config.fish; and rm "$__fish_config_dir/config.fish.bak"; or return 1
    end

    echo "__bdm_fish_functions$marker" >> $__fish_config_dir/config.fish; or return 1
    __bdm_fish_functions; or return 1

    ### Environment Indicator ###
    if not $_bdm_quiet
        set -gx __bdm_major_minor (bdm-config --version | sed -n  '1 s/.*v\([0-9]*.[0-9]*\).*/\1/p')
        if not type -qt __bdm_fish_prompt_original
            functions --copy fish_prompt __bdm_fish_prompt_original
        end
        # NB: overrides prompt for current session only
        function fish_prompt
            echo "[bdm-$__bdm_major_minor] "(__bdm_fish_prompt_original)
        end
    end
end

if source_thisbdm $argv
    _bdm_ok "[OK] You have successfully sourced BioDynaMo's environment."
    _thisbdm_cleanup
else
    _bdm_err "[ERR] BioDynaMo's environment could not be sourced."
    _thisbdm_cleanup
end
