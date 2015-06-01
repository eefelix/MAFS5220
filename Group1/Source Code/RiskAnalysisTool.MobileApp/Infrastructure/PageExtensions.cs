using System;
using System.Collections.Generic;
using System.Text;
using System.Threading.Tasks;
using Xamarin.Forms;

namespace RiskAnalysisTool.MobileApp.Views
{
    public static class PageExtensions
    {
        public static async Task<int?> DisplaySelectionBox(this Page page, string title, params string[] items)
        {
            string selectedItem = await page.DisplayActionSheet(title, "Cancel", null, items);
            if (selectedItem == null || selectedItem == "Cancel")
            {
                return null;
            }
            else
            {
                return Array.IndexOf(items, selectedItem);
            }
        }
    }
}
