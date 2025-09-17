# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: 2020-2025 PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

# coding: utf-8

import base64
import os

import streamlit as st


def use_flexo():
    """
    Apply Flexo font throughout
    """
    font_extensions = [".woff2", ".woff", ".ttf"]

    base_dir = os.path.dirname(os.path.abspath(__file__))
    font_base_path = os.path.join(base_dir, "design", "fonts", "Flexo-Medium")

    for ext in font_extensions:
        font_path = font_base_path + ext
        if os.path.exists(font_path):
            try:
                with open(font_path, "rb") as f:
                    # Since streamlit often serves font files with the wrong MIME type,
                    # which results in decoding problems on fetching later, the solution
                    # here is to convert the font file to base64 data and embed it
                    # directly into the custom CSS
                    font_data = base64.b64encode(f.read()).decode("utf-8")

                # Establish a MIME type to ensure the correct type is used by the browser
                mime_type = {
                    ".woff2": "font/woff2",
                    ".woff": "font/woff",
                    ".ttf": "font/ttf",
                }[ext]

                st.markdown(
                    f"""
                    <style>
                    @font-face {{
                        font-family: 'Flexo';
                        src: url(data:{mime_type};base64,{font_data}) format('{ext[1:]}');
                        font-weight: normal;
                        font-style: normal;
                    }}

                    body, p, div, span, pre, .stMarkdown, .stTitle, .stText {{
                        font-family: 'Flexo', sans-serif;
                    }}

                    h1, h2, h3, h4, h5, h6 {{
                        font-family: 'Flexo', sans-serif !important;
                    }}

                    .plot-container {{
                        font-family: 'Flexo', sans-serif;
                    }}
                    </style>
                    """,
                    unsafe_allow_html=True,
                )

            except Exception as e:
                st.Error(f"Error loading font: {e}")


def apply_sidebar_styles():
    """
    Style the Page navigation and Parameters part of the sidebar
    """
    st.markdown(
        """
        <style>
        div[data-testid="stSidebarNav"]::before {
            content: "Page Navigation";
            margin-left: 20px;
            margin-bottom: 12px;
            font-size: 1.2em;
            font-weight: 600;
            position: relative;
            top: 4px;
            display: block;
        }
        /* Make the navigation divider match st.divider() style (shorter) */
        div[data-testid="stSidebarNavSeparator"] {
            margin: 0.68em 1.2rem 0.5rem;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )


def apply_sidebar_chart_nav_styles():
    """
    Style the Chart navigation part of the sidebar
    """
    st.markdown(
        """
            <p style='font-size: 1.2em; font-weight: 600; margin-bottom: -4px;'>
            Chart navigation
            </p>
            """,
        unsafe_allow_html=True,
    )
    st.markdown(
        """
            <style>
            .nav-link {
                display: flex;
                align-items: center;
                padding-left: 16px;
                text-decoration: none !important;
                color: inherit !important;
                border-radius: 6px;
                text-align: left;
                width: 100%;
                box-sizing: border-box;
                transition: all 0.2s ease !important;
                font-family: inherit !important;
                line-height: 2;
                margin-bottom: -24px;
                background-color: transparent;
            }
            .nav-link:hover {
                background-color: rgba(128, 128, 128, 0.1) !important;
                color: inherit !important;
            }
            /* Dark mode support */
            @media (prefers-color-scheme: dark) {
                .nav-link:hover {
                    background-color: rgba(255, 255, 255, 0.1) !important;
                }
            }
            </style>
            """,
        unsafe_allow_html=True,
    )


def apply_radio_menu_styles():
    st.markdown(
        """
        <style>
        /* Centre radio buttons (for bar_with_filter) */
        div[data-testid="stVerticalBlockBorderWrapper"] div.stRadio > div[role="radiogroup"] {
            display: flex;
            justify-content: center;
            gap: 1rem;
        }
        /* Centre date filter elements (for graphs with date filters) */
        div[data-testid="stElementContainer"] {
            display: flex;
            justify-content: center;
        }
        div.stButtonGroup {
            width: auto;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )
